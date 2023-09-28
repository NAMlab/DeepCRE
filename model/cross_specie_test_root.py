from tensorflow.keras import models
import os
import pandas as pd
from utils import prepare_valid_seqs
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'

mapped_read_counts = ['zea_root_counts.csv', 'solanum_root_counts.csv', 'arabidopsis_root_counts.csv',
                      'sbicolor_root_counts.csv']
gene_models = ['Zea_mays.Zm-B73-REFERENCE-NAM-5.0.52.gtf', 'Solanum_lycopersicum.SL3.0.52.gtf',
               'Arabidopsis_thaliana.TAIR10.52.gtf', 'Sorghum_bicolor.Sorghum_bicolor_NCBIv3.52.gtf']
genomes = ['Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa', 'Solanum_lycopersicum.SL3.0.dna.toplevel.fa',
           'Arabidopsis_thaliana.TAIR10.dna.toplevel.fa', 'Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa']
chrom_nums = {'Arabidopsis': 5, 'Solanum': 12, 'Zea': 10, 'Sorghum': 10}

cross_specie_results = []
for tpm, gene_model, genome in zip(mapped_read_counts, gene_models, genomes):
    chroms = [str(x) for x in range(1, chrom_nums[genome.split('_')[0]]+1)]
    x_val, y_val, _ = prepare_valid_seqs(genome, gene_model, tpm, chroms)
    for model in os.listdir('saved_models'):
        if model.split('_')[0] != tpm.split('_')[0] and model.endswith('terminatorroot.h5'):
            print(model)
            loaded_model = models.load_model(f'saved_models/{model}')
            evaluation = loaded_model.evaluate(x_val, y_val)
            print(genome.split('_')[0], model.split('_')[0])
            cross_specie_results.append([model.split('_')[0], genome.split('_')[0], evaluation[-1]])

cross_specie_results = pd.DataFrame(cross_specie_results,
                                    columns=['train_specie', 'cross_specie', 'accuracy'])
print(cross_specie_results.head())
cross_specie_results.to_csv('../results/cross_specie_root_result.csv', index=False)