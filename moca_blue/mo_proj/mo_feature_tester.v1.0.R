library(tidyverse)
##############################################################
PROJECT <- "_on_SolaITAG4_ch01_0e3-cwm"
SPEC <- "Soly"
MODEL <- "MSR"
DATE <- "20230703"
##############################################################
#file_path <- "GO_enrichment_result_table.csv" 
#file_path1 <-"feat_enrichment_result_table.csv"    #Results
##############################################################
dirpath_1 <- "../ref_seq"
dirpath_2 <- "./out"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
file1 <- "SolyMSR_on_Spe-ch01-0e3_gene_none20230530feat_mima_q1q9.csv"
file3 <- "mapman_sopen.txt"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
file_path_out <- file.path(dirpath_2, paste0(DATE,"_",PROJECT,"_mo-feat"))
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
FILTER<- "q1q9"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
mm0 <- read.table(
  file.path(
    dirpath_2,
    file1),
  header=TRUE,
  sep=",")

mapman <- read.table(file3,
                     header=TRUE,
                     sep="\t", quote = "")

########### ############## ############ #####################
colnames(mapman)[colnames(mapman)=="IDENTIFIER"]<-"loc_ID"

mm0$loc_ID <- tolower(mm0$loc_ID)
mm0$loc_ID <- gsub("[-']+", "", mm0$loc_ID)
mm0$loc_ID <- gsub("\\..*", "", mm0$loc_ID) #CAREFULL WITH THE DOTS

mapman$loc_ID <- tolower(mapman$loc_ID)
mapman$loc_ID <- gsub("[-']+", "", mapman$loc_ID)
mapman$loc_ID <- gsub("\\..*", "", mapman$loc_ID) #CAREFULL WITH THE DOTS

mapman0 <- mapman[, c("loc_ID",
                      "BINCODE",
                      "NAME")]

mm0_mapman <- merge(mm0, mapman0, by= "loc_ID")
# Internal check for the number of rows
if (nrow(mm0_mapman) < 2) {
  error_message <- "Please check identifiers chr#:start-end for fasta and gff input"
  stop(error_message)
}
########### ############## ############ ##########
mm0_mapman$BINCODE <- gsub("'", "", mm0_mapman$BINCODE)
mm0_mapman <- mm0_mapman[mm0_mapman$BINCODE != 35.2, ]
mm0_mapman <- mm0_mapman[mm0_mapman$BINCODE != 35.1, ]
mm0_mapman <- mm0_mapman[mm0_mapman$BINCODE != 35, ]
#mm0_mapman$BINCODE <- gsub("\\.([^']*)$", "", mm0_mapman$BINCODE)
#now all is gone after the first . POINT !!!
########### ############## ############ ##########
col_idx <- mm0_mapman[, c("BINCODE", "NAME")]
col_idx <- unique(col_idx)
##################################################
cont_table_A <- table(mm0_mapman$motif, mm0_mapman$BINCODE)
testA <- as.data.frame(cont_table_A)
tfestA <- testA %>%
  pivot_wider(names_from = Var2, values_from = Freq)
##################################################
row_idx <- as.data.frame(tfestA[, c("Var1")])
row_idx <- unique(row_idx)
idx1 <- nrow(row_idx)
row_idx$number <- rownames(row_idx)
colnames(row_idx)[1] <- "motif"
# Exclude the first column from tfestA
tfestA <- tfestA[, -1]
# Convert factor columns to numeric
tfestA <- as.data.frame(lapply(tfestA, as.numeric))
# Assuming your dataframe is called 'tfestA'
sample_size <- sum(tfestA)
results <- vector("list", nrow(tfestA))

#
#for (i in 1:nrow(tfestA)) {
#  row_result <- vector("list", ncol(tfestA))
#  eval1 <- sum(tfestA[i, ]) * sum(tfestA[i,]) / sample_size
#  calc1 <- (tfestA[i,]-eval1)^2/eval1
#  row_result[[i]] <- calc1 
#  results[[i]] <- calc1
#}
#

results <- vector("list", nrow(tfestA))

for (i in 1:nrow(tfestA)) {
  row_result <- vector("list", ncol(tfestA))
  eval1 <- sum(tfestA[i, ]) * sum(tfestA[, i]) / sample_size
  calc1 <- (tfestA[i,]-eval1)^2/eval1
  row_result[[i]] <- calc1 
  results[[i]] <- calc1
}

print(sample_size)
print(min(tfestA))
print(max(tfestA))
print(eval1)

result_matrix <- as.data.frame(do.call(rbind, results))
result_matrix$number <- c(1:idx1)
result_matrix0 <- merge(row_idx,
                        result_matrix,
                        by = "number")


result_df <- pivot_longer(result_matrix0,
                          cols = -c(number, motif),
                          names_to = "column names",
                          values_to = "values")

result_df <- result_df[, -1]
colnames(result_df)[2] <- "BINCODE"
colnames(result_df)[3] <- "x_square"
result_df <- result_df %>%
  mutate(BINCODE = str_replace(BINCODE, "X", ""))
result_df <- result_df %>%
  mutate(BINCODE = gsub("^\\.+|\\.+\\$", "", BINCODE))
result_df <- result_df %>%
  mutate(BINCODE = gsub("^\\.+|\\.+?$", "", BINCODE))
#########################################################################
sorted_results <- result_df[order(result_df$x_square), ]
#sorted_results <- result_df[result_df$x_square != 35.2, ]
top_100_rows <- sorted_results %>%
  arrange(desc(x_square)) %>%
  head(100)
top_100_rows0 <- merge(top_100_rows, col_idx, by= "BINCODE")


#write.csv(top_100_rows0, file = file_path, row.names = FALSE)

write.csv(top_100_rows0, file=paste0(file_path_out,
                              "-goterm.csv"), row.names=FALSE)
#########################################################################
################################################################## ######
#Replace motifs and GO terms in top 100
#Print - export
#Alternative strategy: observed value - expected value/ expected value, test for independence
################################################################## ######

# # # # # # # # # # #  # # # # # # 

unique_values <- unique(mm0_mapman$type.y)
print(unique_values)

cont_table_B <- table(mm0_mapman$motif, mm0_mapman$type.y)

testB =as.data.frame(cont_table_B)
tfestB = testB %>%
  pivot_wider(names_from = Var2, values_from = Freq)



##########################################################################################


##########################################################################################
##########################################################################################

# Check if "mRNA" and "untranscr" columns are present in the data frame
if ("mRNA" %in% colnames(tfestB) && "untranscr" %in% colnames(tfestB)) {
  tfestB$transcr_class <- apply(tfestB[, c("mRNA", "untranscr")], 1, function(row) {
    result <- chisq.test(row, p = c(0.5, 0.5))
    result$p.value
  })
  tfestB$transcr_pref <- ifelse(tfestB$mRNA < tfestB$untranscr, "untranscribed", "transcribed")
} else {
  print("Error: 'mRNA' and/or 'untranscr' columns are not present in the data frame.")
}
##########################################################################################
# Check if "exon" and "intron" columns are present in the data frame
if ("exon" %in% colnames(tfestB) && "intron" %in% colnames(tfestB)) {
  tfestB$feature_class <- apply(tfestB[, c("exon", "intron")], 1, function(row) {
    result <- chisq.test(row, p = c(0.5, 0.5))
    result$p.value
  })
  tfestB$feature_pref <- ifelse(tfestB$exon < tfestB$intron, "intronic", "exonic")
} else {
  print("Error: 'exon' and/or 'intron' columns are not present in the data frame.")
}
########################################################################################## Comparison might be unfair UTR/CDS better
# Check if "CDS" and "intron" columns are present in the data frame
if ("CDS" %in% colnames(tfestB) && "intron" %in% colnames(tfestB)) {
  tfestB$transl_class <- apply(tfestB[, c("CDS", "intron")], 1, function(row) {
    result <- chisq.test(row, p = c(0.5, 0.5))
    result$p.value
  })
  tfestB$transl_pref <- ifelse(tfestB$CDS < tfestB$intron, "intronic", "codogenic")
} else {
  # Handle the case when "CDS" and/or "intron" columns are not present
  # Print an error message or perform alternative actions
  print("Error: 'CDS' and/or 'intron' columns are not present in the data frame.")
}
##########################################################################################
# Check if "CDS" and "intron" columns are present in the data frame
if ("CDS" %in% colnames(tfestB) && "UTR" %in% colnames(tfestB)) {
  tfestB$transl_class <- apply(tfestB[, c("CDS", "UTR")], 1, function(row) {
    result <- chisq.test(row, p = c(0.5, 0.5))
    result$p.value
  })
  tfestB$transl_pref <- ifelse(tfestB$CDS < tfestB$intron, "UTR", "codogenic")
} else {
  # Handle the case when "CDS" and/or "intron" columns are not present
  # Print an error message or perform alternative actions
  print("Error: 'CDS' and/or 'UTR' columns are not present in the data frame.")
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

##########################################################################################
##########################################################################################

write.csv(tfestB, file=paste0(file_path_out,
                              "-features.csv"), row.names=FALSE)
##########################################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Remove columns with sum less than 5
########### ############## ############ ##########
