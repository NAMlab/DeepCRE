import pandas as pd
from utils import onehot
from pyfaidx import Fasta
from tensorflow.keras import models
from tensorflow.keras import backend
import os
import numpy as np
import deeplift
from deeplift.util import get_shuffle_seq_ref_function
from deeplift.dinuc_shuffle import dinuc_shuffle
from deeplift.conversion import kerasapi_conversion as kc
from sklearn.metrics import accuracy_score
from deeplift.util import get_hypothetical_contribs_func_onehot
import h5py
import modisco
from importlib import reload
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
backend.clear_session()
if not os.path.exists('modisco_multi'):
    os.mkdir('modisco_multi')


def compute_actual_and_hypothetical_scores(x_val, y_val, specie, chrom):
    saved_name = f'saved_models/multi_specie_model_{specie}_{chrom}.h5'
    loaded_model = models.load_model(saved_name)
    predicted_class = loaded_model.predict(x_val)
    predicted_class = predicted_class > 0.5
    print(f'{specie} Accuracy {accuracy_score(y_val, predicted_class)}')

    x = []
    for idx, seq in enumerate(x_val):
        if predicted_class[idx] == 1 and y_val[idx] == 1:
            x.append(seq)
        elif predicted_class[idx] == 0 and y_val[idx] == 0:
            x.append(seq)

    x = np.array(x)

    print(f'Number of correct predictions {x.shape[0]}')
    # ---------- Computing importance and hypothetical scores-------------------------------------------#
    deeplift_model = kc.convert_model_from_saved_files(saved_name,
                                                       nonlinear_mxts_mode=deeplift.layers.NonlinearMxtsMode.DeepLIFT_GenomicsDefault)

    deeplift_contribs_func = deeplift_model.get_target_contribs_func(find_scores_layer_idx=0,
                                                                     target_layer_idx=-2)

    contribs_many_refs_func = get_shuffle_seq_ref_function(
        score_computation_function=deeplift_contribs_func,
        shuffle_func=dinuc_shuffle)

    multipliers_func = deeplift_model.get_target_multipliers_func(find_scores_layer_idx=0,
                                                                  target_layer_idx=-2)
    hypothetical_contribs_func = get_hypothetical_contribs_func_onehot(multipliers_func)

    # Once again, we rely on multiple shuffled references
    hypothetical_contribs_many_refs_func = get_shuffle_seq_ref_function(
        score_computation_function=hypothetical_contribs_func,
        shuffle_func=dinuc_shuffle)

    actual_scores = np.squeeze(np.sum(contribs_many_refs_func(task_idx=0,
                                                              input_data_sequences=x,
                                                              num_refs_per_seq=10,
                                                              batch_size=50,
                                                              progress_update=4000), axis=2))[:, :, None] * x

    hyp_scores = hypothetical_contribs_many_refs_func(task_idx=0,
                                                      input_data_sequences=x,
                                                      num_refs_per_seq=10,
                                                      batch_size=50,
                                                      progress_update=4000)

    return x, actual_scores, hyp_scores


def run_modisco(specie):
    save_file = f"modisco_multi/{specie}_modisco.hdf5"
    os.system(f'rm -rf {save_file}')

    h5_data = h5py.File(f'modisco_multi/{specie}_scores.h5', 'r')
    contribution_scores = h5_data.get('contrib_scores')
    hypothetical_scores = h5_data.get('hypothetical_scores')
    one_hots = h5_data.get('one_hots')

    print('contributions', contribution_scores.shape)
    print('hypothetical contributions', hypothetical_scores.shape)
    print('correct predictions', one_hots.shape)
    # -----------------------Running modisco----------------------------------------------#
    # Uncomment to refresh modules for when tweaking code during development:
    reload(modisco.util)
    reload(modisco.pattern_filterer)
    reload(modisco.aggregator)
    reload(modisco.core)
    reload(modisco.seqlet_embedding.advanced_gapped_kmer)
    reload(modisco.affinitymat.transformers)
    reload(modisco.affinitymat.core)
    reload(modisco.affinitymat)
    reload(modisco.cluster.core)
    reload(modisco.cluster)
    reload(modisco.tfmodisco_workflow.seqlets_to_patterns)
    reload(modisco.tfmodisco_workflow)
    reload(modisco)

    null_per_pos_scores = modisco.coordproducers.LaplaceNullDist(num_to_samp=5000)
    tfmodisco_results = modisco.tfmodisco_workflow.workflow.TfModiscoWorkflow(
        # Slight modifications from the default settings
        sliding_window_size=21,
        flank_size=5,
        target_seqlet_fdr=0.01,
        seqlets_to_patterns_factory=modisco.tfmodisco_workflow.seqlets_to_patterns.TfModiscoSeqletsToPatternsFactory(
            # Note: as of version 0.5.6.0, it's possible to use the results of a motif discovery
            # software like MEME to improve the TF-MoDISco clustering. To use the meme-based
            # initialization, you would specify the initclusterer_factory as shown in the
            # commented-out code below:
            # initclusterer_factory=modisco.clusterinit.memeinit.MemeInitClustererFactory(
            #    meme_command="meme", base_outdir="meme_out",
            #    max_num_seqlets_to_use=10000, nmotifs=10, n_jobs=1),
            trim_to_window_size=10,
            initial_flank_to_add=2,
            final_flank_to_add=0,
            final_min_cluster_size=60,
            # use_pynnd=True can be used for faster nn comp at coarse grained step
            # (it will use pynndescent), but note that pynndescent may crash
            # use_pynnd=True,
            n_cores=40)
    )(
        task_names=['task0'],
        contrib_scores={'task0': contribution_scores},
        hypothetical_contribs={'task0': hypothetical_scores},
        one_hot=one_hots,
        null_per_pos_scores=null_per_pos_scores)

    reload(modisco.util)
    grp = h5py.File(save_file, "w")
    tfmodisco_results.save_hdf5(grp)
    grp.close()


def encode_seq(gene_models, fasta, labels, upstream=1000, downstream=500):
    encoded_seqs, true_labels = [], []
    for chrom, start, end, strand, gene_id, _ in gene_models.values:
        if gene_id in labels.index:
            if strand == '+':
                prom_start, prom_end = start - upstream, start + downstream
                term_start, term_end = end - downstream, end + upstream
                if prom_start > 0 and term_start > 0:
                    seq = str(fasta[chrom][prom_start:prom_end]) + 'N'*20 + str(fasta[chrom][term_start:term_end])

                    if len(seq) == 2 * (upstream + downstream) + 20:
                        encoded_seqs.append(onehot(seq))
                        true_labels.append(labels.loc[gene_id, 'true_target'])

            else:
                prom_start, prom_end = end - downstream, end + upstream
                term_start, term_end = start - upstream, start + downstream
                if prom_start > 0 and term_start > 0:
                    seq = str(fasta[chrom][prom_start:prom_end].reverse.complement) + 'N'*20 +\
                          str(fasta[chrom][term_start:term_end].reverse.complement)

                    if len(seq) == 2 * (upstream + downstream) + 20:
                        encoded_seqs.append(onehot(seq))
                        true_labels.append(labels.loc[gene_id, 'true_target'])

    encoded_seqs = np.array(encoded_seqs)
    true_labels = np.array(true_labels)

    return encoded_seqs, true_labels


gene_model = pd.read_csv('gene_models/zea_sol_ara_sor_52.gtf', sep='\t')
gene_model['specie'] = [x.split('_')[0] for x in gene_model['Chromosome']]
genes_labels = pd.read_csv('tpm_counts/zea_sol_ara_sor.csv', sep='\t', index_col=0)
genome = Fasta('genomes/zea_sol_ara_sor_dna.fa',  as_raw=False, sequence_always_upper=True, read_ahead=10000)
print(gene_model.head())

plant_species = ['Arabidopsis', 'Zea', 'Sorghum', 'Solanum']
for specie in plant_species:
    if not os.path.exists(f"modisco_multi/{specie}_modisco.hdf5"):
        final_results = []
        print(f'Testing on {specie}')
        test_specie = gene_model.copy()
        test_specie = test_specie[test_specie['specie'] == specie]

        actual_scores_all, hypothetical_scores_all, onehot_all = [], [], []
        for chrom_num in test_specie['Chromosome'].unique():
            test_specie_chrom = test_specie.copy()
            test_specie_chrom = test_specie_chrom[test_specie_chrom['Chromosome'] == chrom_num]
            test_seqs, test_targets = encode_seq(test_specie_chrom, genome, genes_labels)

            print(f'Running Deeplift {chrom_num}')
            encoded_x, act_imp, hyp_imp = compute_actual_and_hypothetical_scores(test_seqs, test_targets, specie,
                                                                                 chrom_num)
            actual_scores_all.append(act_imp)
            hypothetical_scores_all.append(hyp_imp)
            onehot_all.append(encoded_x)

        # Save scores in h5 format
        if os.path.isfile(f'modisco_multi/{specie}_scores.h5'):
            os.system(f'rm -rf modisco_multi/{specie}_scores.h5')

        actual_scores_all = np.concatenate(actual_scores_all, axis=0)
        hypothetical_scores_all = np.concatenate(hypothetical_scores_all, axis=0)
        onehot_all = np.concatenate(onehot_all, axis=0)

        h = h5py.File(f'modisco_multi/{specie}_scores.h5', 'w')
        h.create_dataset('contrib_scores', data=actual_scores_all)
        h.create_dataset('hypothetical_scores', data=hypothetical_scores_all)
        h.create_dataset('one_hots', data=onehot_all)
        h.close()

        print('Running modisco')
        run_modisco(specie)
