from tensorflow.keras import models
from tensorflow.keras import backend
import os
import numpy as np
import deeplift
from deeplift.util import get_shuffle_seq_ref_function
from deeplift.dinuc_shuffle import dinuc_shuffle
from deeplift.conversion import kerasapi_conversion as kc
from utils import prepare_valid_seqs
from sklearn.metrics import accuracy_score, roc_auc_score
from deeplift.util import get_hypothetical_contribs_func_onehot
import h5py
import modisco
from importlib import reload
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
backend.clear_session()
if not os.path.exists('modisco'):
    os.mkdir('modisco')


def compute_actual_and_hypothetical_scores(fasta, gtf, tpms, specie):
    actual_scores_all, hypothetical_scores_all, onehot_all = [], [], []
    for saved_model in os.listdir('saved_models'):
        if saved_model.startswith(specie) and saved_model.endswith('terminatorroot.h5'):
            print(saved_model)
            val_chrom = saved_model.split('_')[2]
            x_val, y_val, _ = prepare_valid_seqs(fasta=fasta, gtf=gtf, tpms=tpms, val_chrom=val_chrom)

            saved_name = f'saved_models/{saved_model}'
            loaded_model = models.load_model(saved_name)
            predicted_class = loaded_model.predict(x_val)
            predicted_class = predicted_class > 0.5
            print('Accuracy', accuracy_score(y_val, predicted_class))
            x = []
            for idx, seq in enumerate(x_val):
                if predicted_class[idx] == 0 and y_val[idx] == 0:
                    x.append(seq)
            for idx, seq in enumerate(x_val):
                if predicted_class[idx] == 1 and y_val[idx] == 1:
                    x.append(seq)
            x = np.array(x)

            print(f'Number of correct predictions {x.shape[0]}, {x_val.shape}')
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

            actual_scores_all.append(actual_scores)
            hypothetical_scores_all.append(hyp_scores)
            onehot_all.append(x)

    # Save scores in h5 format
    if os.path.isfile(f'modisco/{specie}_scores_root.h5'):
        os.system(f'rm -rf modisco/{specie}_scores_root.h5')

    actual_scores_all = np.concatenate(actual_scores_all, axis=0)
    hypothetical_scores_all = np.concatenate(hypothetical_scores_all, axis=0)
    onehot_all = np.concatenate(onehot_all, axis=0)

    h = h5py.File(f'modisco/{specie}_scores_root.h5', 'w')
    h.create_dataset('contrib_scores', data=actual_scores_all)
    h.create_dataset('hypothetical_scores', data=hypothetical_scores_all)
    h.create_dataset('one_hots', data=onehot_all)
    h.close()


def run_modisco(specie):
    save_file = f"modisco/{specie}_modisco_root.hdf5"
    os.system(f'rm -rf {save_file}')

    h5_data = h5py.File(f'modisco/{specie}_scores_root.h5', 'r')
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
        max_seqlets_per_metacluster=50000,
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
            n_cores=50)
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


species = ['arabidopsis', 'zea', 'solanum', 'sbicolor']
gene_models = ['Arabidopsis_thaliana.TAIR10.52.gtf', 'Zea_mays.Zm-B73-REFERENCE-NAM-5.0.52.gtf',
               'Solanum_lycopersicum.SL3.0.52.gtf', 'Sorghum_bicolor.Sorghum_bicolor_NCBIv3.52.gtf']
genomes = ['Arabidopsis_thaliana.TAIR10.dna.toplevel.fa', 'Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa',
           'Solanum_lycopersicum.SL3.0.dna.toplevel.fa', 'Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa']
pickle_keys = ['ara', 'zea', 'sol', 'sor']
mapped_read_counts = ['arabidopsis_root_counts.csv', 'zea_root_counts.csv', 'solanum_root_counts.csv',
                      'sbicolor_root_counts.csv']


for plant, fasta_file, gtf_file, pickled_key, counts in zip(species, genomes, gene_models, pickle_keys,
                                                            mapped_read_counts):
    if not os.path.exists(f'modisco/{plant}_modisco_root.hdf5'):
        print(f'Computing contribution and hypothetical contribution scores for {plant}-----------------------------\n')
        compute_actual_and_hypothetical_scores(fasta_file, gtf_file, counts, plant)
        print(f'Running TFMoDisco on {plant}------------------------------------------------------------------------\n')
        run_modisco(plant)
