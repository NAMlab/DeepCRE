import pandas as pd
import numpy as np
import os
from utils import FastaSequenceLoader, ConvNetwork
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
mapped_read_counts = ['zea_counts.csv', 'solanum_counts.csv', 'arabidopsis_counts.csv', 'sbicolor_counts.csv']
gene_models = ['Zea_mays.Zm-B73-REFERENCE-NAM-5.0.52.gtf', 'Solanum_lycopersicum.SL3.0.52.gtf',
               'Arabidopsis_thaliana.TAIR10.52.gtf', 'Sorghum_bicolor.Sorghum_bicolor_NCBIv3.52.gtf']
genomes = ['Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa', 'Solanum_lycopersicum.SL3.0.dna.toplevel.fa',
           'Arabidopsis_thaliana.TAIR10.dna.toplevel.fa', 'Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa']
pickle_keys = ['zea', 'sol', 'ara', 'sor']
num_chromosomes = [10, 12, 5, 10]

if not os.path.isdir('../results'):
    os.mkdir('../results')
if not os.path.isdir('saved_models/size_effect'):
    os.mkdir('saved_models/size_effect')

for m_reads, gene_model, genome, num_chr, p_key in zip(mapped_read_counts, gene_models, genomes, num_chromosomes,
                                                       pickle_keys):
    for inner_flank in [100, 300, 700]:
        print('-----------------------------------------------------------------------------\n')
        print(f"Plant: {m_reads.split('_')[0]} Outer flank: 1000nt inner flank: {inner_flank}nt")
        print('-----------------------------------------------------------------------------\n')

        if not os.path.exists(f"../results/{m_reads.split('_')[0]}_{inner_flank}_utr_result.csv"):
            final_training_output = []
            tpm_counts = pd.read_csv(f'tpm_counts/{m_reads}', index_col=0)
            true_targets = []

            for log_count in tpm_counts['logMaxTPM'].values:
                if log_count <= np.percentile(tpm_counts['logMaxTPM'], 25):
                    true_targets.append(0)
                elif log_count >= np.percentile(tpm_counts['logMaxTPM'], 75):
                    true_targets.append(1)
                else:
                    true_targets.append(2)
            tpm_counts['true_target'] = true_targets
            print(tpm_counts.head())

            for val_chromosome in np.arange(1, num_chr+1):
                fastaloader = FastaSequenceLoader(f'genomes/{genome}', f'gene_models/{gene_model}',
                                                  val_chromosome, pickled_val_ids='validation_genes.pickle',
                                                  pickled_key=p_key, upstream=1000, downstream=inner_flank)
                enc_train, enc_val, train_ids, val_ids = fastaloader.extract_seq()
                convnet = ConvNetwork(enc_train, enc_val, train_ids, val_ids, val_chromosome, tpm_counts,
                                      m_reads.split('_')[0], 'promoter_terminator', inner_flank=inner_flank,
                                      size_effect=True)
                val_acc, val_auroc, specie_name, _, _ = convnet.train_network()
                final_training_output.append([val_acc, val_auroc, specie_name, inner_flank])

            final_training_output = pd.DataFrame(final_training_output, columns=['val_acc', 'val_auROC', 'plant',
                                                                                 'size'])
            final_training_output.to_csv(f"../results/{m_reads.split('_')[0]}_{inner_flank}_utr_result.csv",
                                         index=False)