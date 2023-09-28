import pyranges as pr
import pandas as pd
import numpy as np
from Bio.SeqIO import SeqRecord
from Bio import SeqIO

print('Now generating zea_sol_ara_sor_dna.gtf')
gene_models = ['Zea_mays.Zm-B73-REFERENCE-NAM-5.0.52.gtf', 'Solanum_lycopersicum.SL3.0.52.gtf',
               'Arabidopsis_thaliana.TAIR10.52.gtf', 'Sorghum_bicolor.Sorghum_bicolor_NCBIv3.52.gtf']
num_chromosomes = [10, 12, 5, 10]
gtfs = []
for idx, specie in enumerate(gene_models):
    gene_models = pr.read_gtf(f'gene_models/{specie}', as_df=True)
    gene_models = gene_models[gene_models['Feature'] == 'gene']
    gene_models = gene_models[gene_models['gene_biotype'] == 'protein_coding']
    gene_models = gene_models[['Chromosome', 'Start', 'End', 'Strand', 'gene_id']]
    gene_models = gene_models[gene_models['Chromosome'].isin([str(x) for x in range(1, num_chromosomes[idx]+1)])]
    gene_models['Chromosome'] = specie.split('_')[0] + '_' + gene_models['Chromosome']
    gene_models['specie'] = [x.split('_')[0] for x in gene_models['Chromosome']]
    gtfs.append(gene_models)
    print(gene_models.head())

gtfs = pd.concat(gtfs)
gtfs.to_csv('gene_models/zea_sol_ara_sor_52.gtf', sep='\t', index=False)

# For tpm counts
counts_root = ['zea_root_counts.csv', 'solanum_root_counts.csv', 'arabidopsis_root_counts.csv',
               'sbicolor_root_counts.csv']
counts_leaf = ['zea_counts.csv', 'solanum_counts.csv', 'arabidopsis_counts.csv',
               'sbicolor_counts.csv']

print('Now generating zea_sol_ara_sor_roots.csv and zea_sol_ara_sor.csv for root and leaf respectively')
for tissue, tissue_counts in zip(['leaf', 'root'], [counts_leaf, counts_root]):
    df_counts = []
    for count in tissue_counts:
        counts_df = pd.read_csv(f'tpm_counts/{count}')
        true_targets = []
        print(count)
        print(np.percentile(counts_df['logMaxTPM'], 25))
        print(np.percentile(counts_df['logMaxTPM'], 75))
        for log_count in counts_df['logMaxTPM'].values:
            if log_count <= np.percentile(counts_df['logMaxTPM'], 25):
                true_targets.append(0)
            elif log_count >= np.percentile(counts_df['logMaxTPM'], 75):
                true_targets.append(1)
            else:
                true_targets.append(2)
        counts_df['true_target'] = true_targets
        counts_df = counts_df[counts_df['true_target'].isin([0, 1])]
        df_counts.append(counts_df[['gene_id', 'true_target']])
        print(counts_df.head())

    df_counts = pd.concat(df_counts)
    print(df_counts.head())
    if tissue == 'root':
        df_counts.to_csv('tpm_counts/zea_sol_ara_sor_roots.csv', index=False)
    else:
        df_counts.to_csv('tpm_counts/zea_sol_ara_sor.csv', index=False)

print('Now generating zea_sol_ara_sor_dna.fa')
genomes = ['Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa', 'Solanum_lycopersicum.SL3.0.dna.toplevel.fa',
           'Arabidopsis_thaliana.TAIR10.dna.toplevel.fa', 'Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa']


records = []
for genome in genomes:
    for rec in SeqIO.parse(f'genomes/{genome}', format='fasta'):
        print(f"{genome.split('_')[0]}_"+rec.id)
        print(f"{genome.split('_')[0]}_"+rec.description)
        records.append(SeqRecord(seq=rec.seq, id=f"{genome.split('_')[0]}_"+rec.id,
                                 description=f"{genome.split('_')[0]}_"+rec.description))

SeqIO.write(records, handle='genomes/zea_sol_ara_sor_dna.fa', format='fasta')