library(rtracklayer)
library(tximport)
library(tibble)
library(dplyr)

Tx <- rtracklayer::import('spenn_v2.0_gene_models_annot.gff')
Tx<- as_tibble(Tx)
Tx <- Tx%>%
  filter(type=='mRNA')%>%
  select(Name)
Tx$Parent <- substr(Tx$Name, 1, nchar(Tx$Name)-2)

srr <- list.files('../../RNASEQ/Spenn/')
path <- paste('../../RNASEQ/Spenn/', srr, '/abundance.tsv', sep = "")

kal_out <- tximport(path,
                    type = 'kallisto',
                    countsFromAbundance = "lengthScaledTPM",
                    tx2gene = Tx)

tpm <- kal_out$abundance
colnames(tpm) <- srr
tpm <- as_tibble(tpm, rownames='gene_id')
write.csv(tpm, file = 'Spenn_tpm_counts.csv', row.names = F)

