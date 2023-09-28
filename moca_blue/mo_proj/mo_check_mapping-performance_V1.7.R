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
PROJECT <- "Arth-0e3-cwm-W2q1q9"
SPEC <- "Arth"
MODEL <- "M0"
DATA_ORIGIN <- "HDF5-BLAMM"
##############################################################
##############################################################
dirpath_1 <- "../../ref_seq"
dirpath_2 <- "./out"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
file1 <- "20230823ArthM0Arth-0e3-cwm-W2_gene_none-q1q9.csv"
file2 <- "Arth_M0_predictions.csv" # file with predicted probabilies for expression 0-1
file3 <- "Arabidopsis_thaliana.TAIR10.pep.mercator4.5.txt"  # GO-term enrichment
file4 <- "Arabidopsis_thaliana_TPMs-peleke-etal2023.csv"   # file with measurement of models (TPM, quartile classes 0,1,2)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
file_path_in_file3 <- file.path(dirpath_1, file3)
file_path_in_file4 <- file.path(dirpath_1, file4)
file_path_out <- file.path(dirpath_2, paste0(SPEC,MODEL,PROJECT,"_performance"))
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
library(stats)
library(dplyr)
library(stringr)
library(ggplot2)
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

###################################################
#merg_df01 %>%
#  group_by(motif, loc_ID) %>%
#  summarize(count = n()) %>%
#  group_by(motif) %>%
#  summarize(average_count = mean(count))
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
head(ss_model1)
count_of_high_pred <- sum(model0$prob == 1)
count_of_low_pred <- sum(model0$prob == 0)
count_of_high_expr <- sum(ss_model1$class == "high")
count_of_low_expr <- sum(ss_model1$class == "low")
#head(tfestCDE)
tfestCDF$TPR_p0_expr <- tfestCDF$expr_class0 / (tfestCDF$expr_class0 + tfestCDF$expr_class1)
tfestCDF$TPR_p1_expr <-  tfestCDF$expr_class1 / (tfestCDF$expr_class0 + tfestCDF$expr_class1)
tfestCDF$TPR_p0_prob <- tfestCDF$prob_class0 / (tfestCDF$prob_class0 + tfestCDF$prob_class1)
tfestCDF$TPR_p1_prob <-  tfestCDF$prob_class1 / (tfestCDF$prob_class0 + tfestCDF$prob_class1)
tfestCDF$TPR_TF <- tfestCDF$epm_TPpred / (tfestCDF$epm_TPpred + tfestCDF$epm_TNpred)
#tfestCDF$imp_expr_score <- log2((tfestCDF$expr_class0/count_of_high_expr)/(tfestCDF$expr_class1/count_of_low_expr))
tfestCDF$imp_expr_score <-ifelse(grepl("_p0m", tfestCDF$Var1),  log2((tfestCDF$expr_class0/count_of_high_expr)/(tfestCDF$expr_class1/count_of_low_expr)),
                                                      log2((tfestCDF$expr_class1/count_of_low_expr)/(tfestCDF$expr_class0/count_of_high_expr))
)
#tfestCDF$imp_pred_score <- log2((tfestCDF$expr_class0/count_of_high_pred)/(tfestCDF$expr_class1/count_of_low_pred))
tfestCDF$imp_pred_score <-ifelse(grepl("_p0m", tfestCDF$Var1),  log2((tfestCDF$prob_class0/count_of_high_pred)/(tfestCDF$prob_class1/count_of_low_pred)),
                                                      log2((tfestCDF$prob_class1/count_of_low_pred)/(tfestCDF$prob_class0/count_of_high_pred))
)
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
print(summary_table)
summary_table$mean_occ_in_expression_class <- ifelse(summary_table$mo_metacluster==1, 1-(summary_table$mean_occ_in_expression_class), summary_table$mean_occ_in_expression_class)
summary_table$mean_occ_in_predicted_class <- ifelse(summary_table$mo_metacluster==1, 1-(summary_table$mean_occ_in_predicted_class), summary_table$mean_occ_in_predicted_class)
sorted_tfestCDF <- tfestCDF %>%
  arrange(desc(TPR_TF)) %>%
  head(n = 10)
print(file1)
print(summary_table)
print(sorted_tfestCDF)
##################################################################################################
write.csv(tfestCDF, file=paste0(file_path_out,
                                ".csv"), row.names=FALSE)

##################################################################################################
