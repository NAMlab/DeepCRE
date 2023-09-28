##############################################################
#This script is designed to find motif with significantly 
# better or worse predictability
# go for pairwise comparison using jackard similarity
# go for pool-wise comparison using dice similarity
# use random forest and estimate importance scores
# then check if importance scores behave similar across different species
##############################################################
#setwd("/home/ibg-4/Desktop/Rhome/solanum_motifs")
##############################################################
PROJECT <- "ArthS0-ArthM0_0e3-cwm_mima"
SPEC <- "Arth"
MODEL <- "M0"
DATA_ORIGIN <- "HDF5-BLAMM"
##############################################################
##############################################################
dirpath_1 <- "../../ref_seq"
dirpath_2 <- "./out"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
file1 <- "ArthS0ArthS0-Arth_0e3_gene_none20230811-mima.csv"
file2 <- "Arth_M0_predictions.csv" # file with predicted probabilies for expression 0-1
file3 <- "Arabidopsis_thaliana.TAIR10.pep.mercator4.5.txt"  # GO-term enrichment
file4 <- "Arabidopsis_thaliana_TPMs-peleke-etal2023.csv"   # file with measurement of models (TPM, quartile classes 0,1,2)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
file_path_in_file3 <- file.path(dirpath_1, file3)
file_path_in_file4 <- file.path(dirpath_1, file4)
file_path_out <- file.path(dirpath_2, paste0(PROJECT,"_performance"))
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
library(stats)
library(dplyr)
library(stringr)
library(ggplot2)
#install.packages("tidyverse")
library(tidyverse)
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
mm0 <- read.table(
  file.path(
    dirpath_2,
    file1),
  header=TRUE,
  sep=",")
mapman <- read.table(file_path_in_file3,
                     header=TRUE,
                     sep="\t", quote = "")
model1 <- read.table(
  file_path_in_file4,
  header=TRUE,
  sep=",")
########### ############## ############ #####################
########### ############## ############ #####################
if (!file.exists(file2)) {
  mm0 <- read.table(file.path(dirpath_2, file1), header=TRUE, sep=",")
  unique_loc_ID <- unique(mm0$loc_ID)  # Remove duplicates
  model0 <- data.frame(loc_ID = unique_loc_ID, prob = sample(c(0, 1), size = length(unique_loc_ID), replace = TRUE))
  print("file2 missing")
  file2_state <- c("FALSE")
} else {
  model0 <- read.table(file2, header=TRUE, sep=",")
  colnames(model0) <- c("loc_ID", "prob")
  print("file2 exists")
  file2_state <- c("TRUE")
}
########### ############# ############# ######################
#############################################################
model0 <- model0[, c(1, 2)]
# Step 1: Subset DF2 to include only the genes with probability > 0.8
#high_prob_genes <- model0[model0$prob > 0.01, ]
# Step 2: Merge DF1 and high_prob_genes based on loc_ID
merged_df <- merge(mm0, model0, 
                   by = "loc_ID",
                   all.x = TRUE,
                   all.y = TRUE)
merged_df <- na.omit(merged_df)
merg_df00 <- merged_df[, c("loc_ID", "dist_transc_border", "region", "epm", "gen_mstart", "prob")] # may include "loc"
#subset_df <- merg_df00[c("motif", "prob")]
#merg_df00 <-merg_df00[!duplicated(merg_df00), ]
################################################################################ HOW TO SCREW DOWN FINE WORK OF ANNOTATION ...
# FIND A WORKING ALTERNATIVE BY THE TIME !!!!
# extract the required string from IDENTIFIER
colnames(mapman)[colnames(mapman)=="IDENTIFIER"]<-"loc_ID"
merg_df00$loc_ID <- tolower(merg_df00$loc_ID)
mapman$loc_ID <- tolower(mapman$loc_ID)
mapman$loc_ID <- gsub("[-']+", "", mapman$loc_ID)
mapman$loc_ID <- gsub("\\..*", "", mapman$loc_ID) #CAREFULL WITH THE DOTS
mapman$loc_ID <- gsub("[-']+", "", mapman$loc_ID)
mapman$loc_ID <- gsub("_.*", "", mapman$loc_ID) # ... AND THE UNDERSCORES
merg_df00$loc_ID <- gsub("[-']+", "", merg_df00$loc_ID)
merg_df00$loc_ID <- gsub("\\..*", "", merg_df00$loc_ID) #CAREFULL WITH THE DOTS
merg_df00map <- merge(
  merg_df00, mapman,
  by = c("loc_ID"),
  all = FALSE, ignore.case = TRUE)
# remove rows with NA cells
################################################################################
ss_model1 <- model1[, c("gene_id", "logMaxTPM", "true_target")]
ss_model1 <- subset(ss_model1, true_target !=  2)                                 # !!! CHANGES DATASET SIZE BY HALF !!! # Were not included within the model training.set
#ss_model1 <- subset(ss_model1, true_target %in% c(1, 2)) # GENES ARE CLASSIFIED AS 0 (low expression), 1 (HIGH) AND 2 (MEDIUM).
# To EVALUATE the predictiveness of EPMs, genes flagged with 0 are handled as LOW expression genes. 
ss_model1$gene_id <- tolower(ss_model1$gene_id)
ss_model1$gene_id <- gsub("[-']+", "", ss_model1$gene_id)
ss_model1$gene_id <- gsub("\\..*", "", ss_model1$gene_id)
colnames(ss_model1) <- c("loc_ID","logMaxTPM","class")
ss_model1$class <- ifelse(
  ss_model1$class == 0,
  "low", "high")

merg_df01mm<- merge(
  merg_df00map, ss_model1,
  by = c("loc_ID"),
  all = FALSE, ignore.case = TRUE)

merg_df01mm$prob_class <- ifelse(
  merg_df01mm$prob <= 0.5,
  "low", "high")

merg_df01mm$pred_perf<- ifelse(
  merg_df01mm$class == merg_df01mm$prob_class,
  "TRUE", "FALSE")

merg_df01mm$epm_pred_perf <- ifelse(
merg_df01mm$pred_perf == FALSE, "NA",
ifelse(
  grepl("p0m", merg_df01mm$epm) & merg_df01mm$class== "high", "TRUE_high",
  ifelse(
    grepl("p1m", merg_df01mm$epm) & merg_df01mm$class== "low", "TRUE_low", "NA")
)
)
################################################################################
a0_fd <- merg_df01mm[, c(5,1,2,3,4,11,6,12,13,14)]
# Contingency table of "epm" and "loc_ID"
#contingency_table_loc_ID <- table(a0_fd$epm, a0_fd$loc_ID)
# Contingency table of "epm" and "GO"
#contingency_table_GO <- table(merg_df01mm$epm, merg_df01mm$NAME)
# Contingency table of "epm" and "class"
contingency_table_expr_class <- table(a0_fd$epm, a0_fd$class)
# Contingency table of "epm" and "prob"
contingency_table_prob_class <- table(a0_fd$epm, a0_fd$prob_class)  # no value if file2 is missing
# Contingency table of "epm" and "performance"
#contingency_table_pred_perf <- table(a0_fd$epm, a0_fd$pred_perf)    # no value if file2 is missing
# Contingency table of "epm" and "performance"
a0_fd2 <- subset(merg_df01mm, epm_pred_perf!= "NA")
contingency_table_epm_pred_perf <- table(a0_fd2$epm, a0_fd2$epm_pred_perf)    # no value if file2 is missing
################################################################################
################################################################################
#testA =as.data.frame(contingency_table_loc_ID)
#tfestA = testA %>%
#  pivot_wider(names_from = Var2, values_from = Freq)
#testB =as.data.frame(contingency_table_GO)
#tfestB = testB %>%
#  pivot_wider(names_from = Var2, values_from = Freq) # transpose A & B, sum, remove greater than 1, 
################################################################################
#tra_tfestA <- as.data.frame(t(tfestA))
#colnames(tra_tfestA) <- tra_tfestA[1, ]
#tra_tfestA <- tra_tfestA[-1, ]
#colnames(tra_tfestA) <- sub("^.*?(Soly_M\\d+_p\\d+m\\d+F).*", "\\1", colnames(tra_tfestA))
#colnames(tra_tfestA) <- sub("^.*?(Soly_M\\d+_p\\d+m\\d+R).*", "\\1", colnames(tra_tfestA))
#tra_tfestA[, -1] <- sapply(tra_tfestA[, -1], as.numeric)
#tra_tfestA[tra_tfestA > 0] <- TRUE
#tra_tfestA0 <- tra_tfestA[rowSums(tra_tfestA == 1) >= 3, ]
#tra_tfestA0[tra_tfestA0 > 0] <- 1
#filtered_df <- tra_tfestA
#filtered_df[filtered_df > 0] <- 1
#filtered_df <- sapply(filtered_df, as.numeric)
#sums <- colSums(filtered_df, na.rm = TRUE)
#sum_df <- data.frame(columns = colnames(filtered_df), sums = sums)

#plot <- ggplot(sum_df, aes(x = columns, y = sums, color = columns)) +
#  geom_point(size = 3) +
#  labs(title = "Column sums", x = "epms", y = "sums") +
#  theme_minimal() +
#  theme(legend.position = "bottom")

# Save the plot as a PDF file
#ggsave(file=paste0(file_path_out,
#                  "sum-q1q9.pdf"), plot, width = 8, height = 6)
#######################################################################################
#sums <- colSums(filtered_df, na.rm = TRUE)
#averages <- sums / nrow(filtered_df)

#sum_df <- data.frame(columns = colnames(filtered_df), sums = sums, averages = averages)
# Create the barplot using ggplot2
#plot2 <- ggplot(sum_df, aes(x = columns, y = averages, color = columns)) +
#  geom_point(size = 3) +
#  labs(title = "Column Averages", x = "Columns", y = "Average") +
#  theme_minimal() +
#  theme(legend.position = "bottom")

# Save the plot as a PDF file
#ggsave(file=paste0(file_path_out,
#                   "average-q1q9.pdf"), plot2, width = 8, height = 6)
#######################################################################################
#head(tra_tfestA0, 10)
#######################################################################################
#######################################################################################
# tra_tfestA0 CONTAINS A RAW SET OF ALL epm MACTHES ACROSS THE DIFFERENT GENES AND CAN BE USED TO FIND BEST COMBINATIONS
#######################################################################################
#######################################################################################
#######################################################################################
#tra_tfestA0$epms <- apply(tra_tfestA0, 1, function(row) {
#  names(row)[row == 1] %>% paste(collapse = "-")
#})
#fmo_comb <- data.frame(row.names = row.names(tra_tfestA0), epms = tra_tfestA0$epms)
#fmo_comb$loc_ID <- row.names(fmo_comb)
#contingency_table_epms <- table(fmo_comb$epms, fmo_comb$loc_ID)
#testF =as.data.frame(contingency_table_epms)
#tfestF = testF %>%
#  pivot_wider(names_from = Var2, values_from = Freq)
#######################################################################################
################################################################################
testC =as.data.frame(contingency_table_expr_class)
tfestC = testC %>%
  pivot_wider(names_from = Var2, values_from = Freq)
colnames(tfestC)<- c("Var1", "expr_class0", "expr_class1")

testD =as.data.frame(contingency_table_prob_class)
tfestD = testD %>%
  pivot_wider(names_from = Var2, values_from = Freq)
colnames(tfestD)<- c("Var1", "prob_class0", "prob_class1")

#testE =as.data.frame(contingency_table_pred_perf)
#tfestE = testE %>%
#  pivot_wider(names_from = Var2, values_from = Freq)

testF =as.data.frame(contingency_table_epm_pred_perf)
tfestF = testF %>%
  pivot_wider(names_from = Var2, values_from = Freq)
################################################################################
tfestCDF <- merge( tfestC, 
                   merge(tfestD, tfestF, by = "Var1"), by = "Var1")

tfestCDF$epm_TPpred <- ifelse(tfestCDF$TRUE_high== 0, tfestCDF$TRUE_low, tfestCDF$TRUE_high)
tfestCDF$epm_TNpred <- c((tfestCDF$expr_class0+tfestCDF$expr_class1)-tfestCDF$epm_TPpred)

tfestCDF <- tfestCDF[, -c(6, 7)] 

tfestCDF$mo_metacluster <- 0
tfestCDF$mo_metacluster[grep("_p1m", tfestCDF$Var1)] <- 1


tfestCDF$prob_rat01 <- ifelse(
  ifelse(tfestCDF$mo_metacluster == 0, 
                              tfestCDF$prob_class0/tfestCDF$prob_class1-1, 
                              tfestCDF$prob_class1/tfestCDF$prob_class0-1) > 0, 1, 0)

tfestCDF$expr_rat01 <- ifelse(
  ifelse(tfestCDF$mo_metacluster == 0, 
                              tfestCDF$expr_class0/tfestCDF$expr_class1-1, 
                              tfestCDF$expr_class1/tfestCDF$expr_class0-1) > 0, 1, 0)
#############################
#### if epm is p1 or p0 and represented in a certain expression class should be calculated here
#head(tfestCDE)
tfestCDF$TPR_p0_expr <- tfestCDF$expr_class0 / (tfestCDF$expr_class0 + tfestCDF$expr_class1)
tfestCDF$TPR_p1_expr <-  tfestCDF$expr_class1 / (tfestCDF$expr_class0 + tfestCDF$expr_class1)
tfestCDF$TPR_p0_prob <- tfestCDF$prob_class0 / (tfestCDF$prob_class0 + tfestCDF$prob_class1)
tfestCDF$TPR_p1_prob <-  tfestCDF$prob_class1 / (tfestCDF$prob_class0 + tfestCDF$prob_class1)
tfestCDF$TPR_TF <- tfestCDF$epm_TPpred / (tfestCDF$epm_TPpred + tfestCDF$epm_TNpred)
#tfestCDE$TNR <- ifelse(tfestCDE$mo_metacluster == 0, tfestCDE$expr_class1 / (tfestCDE$expr_class1 + tfestCDE$expr_class0),
#                       tfestCDE$expr_class0 / (tfestCDE$expr_class1 + tfestCDE$expr_class0))
# !# !# !# THIS TPR/TNR DOES NOT WORK AS INTENDED. # # # # # # # # #  RESULTS LOOK NOT OK ! ! ! ! ! ! ! ! !
#tfestCDE$TPR <- c(tfestCDE$`TRUE`/(tfestCDE$`FALSE`+tfestCDE$`TRUE`))
#tfestCDE$TNR <- c(tfestCDE$`FALSE` / (tfestCDE$`FALSE` + tfestCDE$`TRUE`))
#############################
tfestCDF$chi_expr_class <- apply(tfestCDF[, c("expr_class0", "expr_class1")], 1, function(row) {
  result <- chisq.test(row, p = c(0.5, 0.5))
  result$p.value
})
tfestCDF$chi_prob_class <- apply(tfestCDF[, c("prob_class0", "prob_class1")], 1, function(row) {
  result <- chisq.test(row, p = c(0.5, 0.5))
  result$p.value
})
tfestCDF$chi_epm_TRUE <- apply(tfestCDF[, c("epm_TNpred", "epm_TPpred")], 1, function(row) {
  result <- chisq.test(row, p = c(0.5, 0.5))
  result$p.value
})
###########################################################################
# Aplly fisher exact test?
#############################
if (file2_state == "FALSE") {
  cols_to_replace <- c("prob_class0", "prob_class1", "FALSE", "TRUE", "prob_rat01", "chi_prob_class", "chi_TF")
  tfestCDF[cols_to_replace] <- NA
}
##################################################################################################
summary_table <- tfestCDF %>%
  group_by(mo_metacluster) %>%
  summarize(
    matches = sum(expr_class0 + expr_class1),
    mean_occ_in_expression_class = mean(TPR_p0_expr),
    mean_occ_in_predicted_class = mean(TPR_p0_prob),
    mean_occ_in_true_positives = mean(TPR_TF),
    mean_p.val_occ_in_expression_class = mean(chi_expr_class),
    mean_p.val_occ_in_predicted_class = mean(chi_prob_class),
    mean_p.val_occ_in_true_positives = mean(chi_epm_TRUE)
  )

sorted_tfestCDF <- tfestCDF %>%
  arrange(desc(TPR_TF)) %>%
  head(n = 10)

print(file1)
print(summary_table)
print(sorted_tfestCDF)
##################################################################################################


write.csv(tfestCDF, file=paste0(file_path_out,
                                ".csv"), row.names=FALSE)


