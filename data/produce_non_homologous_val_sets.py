import os
import pickle
import pandas as pd
from Bio import SeqIO
proteins = ['Arabidopsis_thaliana.TAIR10.pep.all.fa', 'Zea_mays.Zm-B73-REFERENCE-NAM-5.0.pep.all.fa',
            'Solanum_lycopersicum.SL3.0.pep.all.fa', 'Sorghum_bicolor.Sorghum_bicolor_NCBIv3.pep.all.fa']

blast_outputs = ['BLAST_ara_to_ara', 'BLAST_zea_to_zea', 'BLAST_sol_to_sol', 'BLAST_sor_to_sor']
validation_genes = dict()
for proteome, blast_output in zip(proteins, blast_outputs):
    if os.path.exists(f'proteomes/{blast_output}'):
        info = []
        for rec in SeqIO.parse(f'proteomes/{proteome}', 'fasta'):
            description = rec.description.split(' ')
            protein_id = description[0]
            gene_id = description[3].split(':')[-1]
            chrom = description[2].split(':')[2]
            info.append([gene_id, protein_id, chrom])
        info = pd.DataFrame(info, columns=['gene_id', 'protein_id', 'chrom'])
        info.index = info.protein_id.tolist()
        print(info.head())
        blast_out = pd.read_csv(f'proteomes/{blast_output}', sep='\t',
                                names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstat', 'qend',
                                       'sstart', 'send', 'evalue', 'bitscore'])
        blast_out['qgene'] = info.loc[blast_out.qseqid.tolist(), 'gene_id'].values
        blast_out['sgene'] = info.loc[blast_out.sseqid.tolist(), 'gene_id'].values
        blast_out['qchrom'] = info.loc[blast_out.qseqid.tolist(), 'chrom'].values
        blast_out['schrom'] = info.loc[blast_out.sseqid.tolist(), 'chrom'].values
        blast_out = blast_out[blast_out['evalue'] < 0.001]
        blast_out = blast_out[blast_out['bitscore'] >= 50]
        blast_out = blast_out[~blast_out['schrom'].isin(['Pt', 'Mt'])]
        blast_out = blast_out[~blast_out['qchrom'].isin(['Pt', 'Mt'])]
        blast_out['homologs_pairs'] = [x + '_' + y for x, y in zip(blast_out.qchrom, blast_out.schrom)]
        val_set = []
        for gene_grp in blast_out.groupby('qgene'):
            if len(gene_grp[1].homologs_pairs.unique()) == 1:
                val_set.append(gene_grp[0])
        validation_genes[f"{blast_output.split('_')[-1]}"] = val_set

if os.path.exists('../model/validation_genes.pickle'):
    os.remove('../model/validation_genes.pickle')
with open('../model/validation_genes.pickle', 'wb') as pickle_file:
    pickle.dump(validation_genes, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)
