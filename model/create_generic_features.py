import pandas as pd
import pyranges as pr
from pyfaidx import Fasta
pd.options.display.width=0


def gc(seq):
    seq = str(seq)
    gc_count = seq.count('G') + seq.count('C')
    return gc_count/len(seq)


def cpg_perc(seq):
    cpg_counts = seq.count('CG')
    return cpg_counts/len(seq) * 100


def get_proximal_prom_and_term_features(gtf, fasta, chroms, flank=1000):
    gene_models = pr.read_gtf(f'gene_models/{gtf}', as_df=True)
    gene_models = gene_models[gene_models['Chromosome'].isin(chroms)]
    fasta = Fasta(f'genomes/{fasta}', as_raw=False, sequence_always_upper=True, read_ahead=10000)
    gene_models = gene_models[gene_models['gene_biotype'] == 'protein_coding']
    gene_models = gene_models[gene_models['Feature'] == 'gene']
    gene_models = gene_models[['Chromosome', 'Start', 'End', 'Strand', 'gene_id']]
    chrom_num, gene_id, prom_gc, term_gc, prom_cpg, term_cpg = [], [], [], [], [], []
    for chrom, start, end, strand, gene in gene_models.values:
        if strand == '-':
            prom_start, prom_end = end, end + flank
            term_start, term_end = start - flank, start
            term_start = 0 if term_start < 0 else term_start
        else:
            prom_start, prom_end = start - flank, start
            prom_start = 0 if prom_start < 0 else prom_start
            term_start, term_end = end, end + flank

        promoter = fasta[chrom][prom_start:prom_end]
        terminator = fasta[chrom][term_start:term_end]

        if strand == '-':
            promoter = promoter.reverse.complement.seq
            terminator = terminator.reverse.complement.seq
        else:
            promoter = promoter.seq
            terminator = terminator.seq
        gene_id.append(gene)
        chrom_num.append(chrom)
        prom_gc.append(gc(promoter))
        term_gc.append(gc(terminator))
        prom_cpg.append(cpg_perc(promoter))
        term_cpg.append(cpg_perc(terminator))

    return gene_id, chrom_num, prom_gc, prom_cpg, term_gc, term_cpg


def get_utr_features(gtf, fasta, chroms):
    utr_length_5, gc_content_5, cpg_5, utr_length_3, gc_content_3, cpg_3, gene_ids = [], [], [], [], [], [], []
    fasta = Fasta(f'genomes/{fasta}', as_raw=False, sequence_always_upper=True, read_ahead=10000)
    gene_models = pr.read_gtf(f'gene_models/{gtf}', as_df=True)
    gene_models = gene_models[gene_models['Chromosome'].isin(chroms)]
    gene_models = gene_models[gene_models['gene_biotype'] == 'protein_coding']
    gene_models = gene_models[['Chromosome', 'Feature', 'Start', 'End', 'Strand', 'gene_id']]
    gene_models = gene_models[gene_models['Feature'].isin(['five_prime_utr', 'three_prime_utr'])]
    for gene in gene_models['gene_id'].unique():
        gene_model_gene = gene_models.copy()
        gene_model_gene = gene_model_gene[gene_model_gene['gene_id'] == gene]
        gene_ids.append(gene)

        for utr_df in gene_model_gene.groupby('Feature'):
            # If no UTR is annotated for the gene
            # For genes with no UTR annotation we consider the UTR to be absent
            if utr_df[1].shape[0] == 0:
                if utr_df[0] == 'five_prime_utr':
                    utr_length_5.append(0)
                    gc_content_5.append(0)
                    cpg_5.append(0)
                elif utr_df[0] == 'three_prime_utr':
                    utr_length_3.append(0)
                    gc_content_3.append(0)
                    cpg_3.append(0)

            # If gene has a UTR annotated
            # some genes have more than one of a specific UTR annotation, so we take the longest
            else:
                utr_df_copy = utr_df[1]
                utr_df_copy['utr_lengths'] = utr_df_copy['End'] - utr_df_copy['Start']
                utr_df_copy = utr_df_copy[utr_df_copy['utr_lengths'] == utr_df_copy['utr_lengths'].max()]
                chrom = utr_df_copy.Chromosome.values[0]
                utr_start = utr_df_copy.Start.values[0]
                utr_end = utr_df_copy.End.values[0]
                strand = utr_df_copy.Strand.values[0]
                utr_length = utr_end - utr_start
                utr_sequence = fasta[chrom][utr_start:utr_end]
                if strand == '-':
                    utr_sequence = utr_sequence.reverse.complement.seq
                else:
                    utr_sequence = utr_sequence.seq
                if utr_df[0] == 'five_prime_utr':
                    utr_length_5.append(utr_length)
                    gc_content_5.append(gc(utr_sequence))
                    cpg_5.append(cpg_perc(utr_sequence))
                elif utr_df[0] == 'three_prime_utr':
                    utr_length_3.append(utr_length)
                    gc_content_3.append(gc(utr_sequence))
                    cpg_3.append(cpg_perc(utr_sequence))
    return gene_ids, utr_length_5, gc_content_5, cpg_5, utr_length_3, gc_content_3, cpg_3


gtfs = ['Zea_mays.Zm-B73-REFERENCE-NAM-5.0.52.gtf', 'Solanum_lycopersicum.SL3.0.52.gtf',
        'Arabidopsis_thaliana.TAIR10.52.gtf', 'Sorghum_bicolor.Sorghum_bicolor_NCBIv3.52.gtf']
genomes = ['Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa', 'Solanum_lycopersicum.SL3.0.dna.toplevel.fa',
           'Arabidopsis_thaliana.TAIR10.dna.toplevel.fa', 'Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa']
num_chromosomes = [10, 12, 5, 10]

for gtf_file, genome, num_chroms in zip(gtfs, genomes, num_chromosomes):
    chromosomes = [str(x) for x in range(1, num_chroms+1)]
    print(gtf_file)

    genes, utr_5, gc_5, cpg5, utr_3, gc_3, cpg3 = get_utr_features(gtf_file, genome, chromosomes)
    df_utr_feats = pd.DataFrame({
        'gene_id': genes,
        "5'UTR length": utr_5,
        "GC 5'UTR ": gc_5,
        "CpG 5'UTR": cpg5,
        "3'UTR length": utr_3,
        "GC 3'UTR ": gc_3,
        "CpG 3'UTR": cpg3,
            })
    print(df_utr_feats.head())

    g_ids, chrom_nums, p_gc, p_cpg, t_gc, t_cpg = get_proximal_prom_and_term_features(gtf_file, genome, chromosomes)
    df_prox_prom_term = pd.DataFrame({
        'gene_id': g_ids,
        'Chromosome': chrom_nums,
        'GC promoter': p_gc,
        'CpG promoter': p_cpg,
        'GC terminator': t_gc,
        'CpG terminator': t_cpg
    })
    print(df_prox_prom_term.head())

    generated_features = df_prox_prom_term.merge(df_utr_feats, how='inner', on='gene_id')
    generated_features.to_csv(f"{gtf_file.split('_')[0]}_generated_features.csv")
