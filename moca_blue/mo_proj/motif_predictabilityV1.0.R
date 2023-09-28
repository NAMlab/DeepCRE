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
PROJECT <- "Sly_mo_on_Spe_0e3_percentiles"
SPEC <- "Spenn"
MODEL <- "MSR"
DATA_ORIGIN <- "motif_matches"
DATE <- "20230530"
##############################################################
##############################################################
dirpath_1 <- "../ref_seq"
dirpath_2 <- "./out"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
file1 <- "SolyMSR_on_Spe-ch01-0e3_gene_none20230530-mima.csv"
file2 <- "msr_predictions_on_pennellii.csv"
file3 <- "mapman_sopen.txt"
file4 <- "spennellii.csv"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
file_path_out <- file.path(dirpath_2, paste0(DATE,"_",PROJECT,"_mo-predict-map"))
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
library(stats)
library(dplyr)
library(stringr)
library(ggplot2)
#library(ape)
#library(hrbrthemes)
#library(randomForestSRC)
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
mm0 <- read.table(
  file.path(
    dirpath_2,
    file1),
  header=TRUE,
  sep=",")
model0 <- read.table(file2,
  header=TRUE,
  sep=",")
mapman <- read.table(file3,
                     header=TRUE,
                     sep="\t", quote = "")
model1 <- read.table(
  file4,
  header=TRUE,
  sep=",")
colnames(model0) <- c("loc_ID","prob")
########### ############## ############ #####################
########### ############## ############ #####################
model0 <- model0[, c(1, 2)]
# Step 1: Subset DF2 to include only the genes with probability > 0.8
#high_prob_genes <- model0[model0$prob > 0.01, ]
# Step 2: Merge DF1 and high_prob_genes based on loc_ID
merged_df <- merge(mm0, model0, 
                   by = "loc_ID",
                   all.x = TRUE,
                   all.y = TRUE)
merged_df <- na.omit(merged_df)
merg_df00 <- merged_df[, c(1, 6, 7, 8,9, 13)]
subset_df <- merg_df00[c("motif", "prob")]
################################################################################
# extract the required string from IDENTIFIER
colnames(mapman)[colnames(mapman)=="IDENTIFIER"]<-"loc_ID"
merg_df00$loc_ID <- tolower(merg_df00$loc_ID)
mapman$loc_ID <- tolower(mapman$loc_ID)
mapman$loc_ID <- gsub("[-']+", "", mapman$loc_ID)
mapman$loc_ID <- gsub("\\..*", "", mapman$loc_ID) #CAREFULL WITH THE DOTS
merg_df00map <- merge(
  merg_df00, mapman,
  by = c("loc_ID"),
  all = FALSE, ignore.case = TRUE)
# remove rows with NA cells
################################################################################
ss_model1 <- model1[, c("gene_id", "logMaxTPM", "true_target")]
ss_model1 <- subset(ss_model1, logMaxTPM !=  0)
#ss_model1 <- subset(ss_model1, true_target %in% c(1, 2)) # GENES ARE CLASSIFIED AS 0 (NO EXPR), 1 (LOW <0.5) AND HIGH (...).
# To EVALUATE the predictiveness of EPMs, genes falgged with 0 are handled as LOW expression genes. 
ss_model1$gene_id <- tolower(ss_model1$gene_id)
ss_model1$gene_id <- gsub("[-']+", "", ss_model1$gene_id)
ss_model1$gene_id <- gsub("\\..*\\.", ".", ss_model1$gene_id)
colnames(ss_model1) <- c("loc_ID","logMaxTPM","class")
ss_model1$class <- ifelse(
  ss_model1$class == 2,
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
################################################################################
a0_fd <- merg_df01mm[, c(5,1,2,3,4,11,6,12,13,14)]
# Contingency table of "motif" and "loc_ID"
contingency_table_loc_ID <- table(a0_fd$motif, a0_fd$loc_ID)
# Contingency table of "motif" and "GO"
contingency_table_GO <- table(merg_df01mm$motif, merg_df01mm$NAME)
# Contingency table of "motif" and "class"
contingency_table_expr_class <- table(a0_fd$motif, a0_fd$class)
# Contingency table of "motif" and "prob"
contingency_table_prob_class <- table(a0_fd$motif, a0_fd$prob_class)
# Contingency table of "motif" and "performance"
contingency_table_pred_perf <- table(a0_fd$motif, a0_fd$pred_perf)
################################################################################
################################################################################
testA =as.data.frame(contingency_table_loc_ID)
tfestA = testA %>%
  pivot_wider(names_from = Var2, values_from = Freq)
testB =as.data.frame(contingency_table_GO)
tfestB = testB %>%
  pivot_wider(names_from = Var2, values_from = Freq) # transpose A & B, sum, remove greater than 1, 
################################################################################

tra_tfestA <- as.data.frame(t(tfestA))
colnames(tra_tfestA) <- tra_tfestA[1, ]
tra_tfestA <- tra_tfestA[-1, ]
colnames(tra_tfestA) <- sub("^.*?(Soly_M\\d+_p\\d+m\\d+F).*", "\\1", colnames(tra_tfestA))
colnames(tra_tfestA) <- sub("^.*?(Soly_M\\d+_p\\d+m\\d+R).*", "\\1", colnames(tra_tfestA))
tra_tfestA[, -1] <- sapply(tra_tfestA[, -1], as.numeric)
#tra_tfestA[tra_tfestA > 0] <- TRUE
#tra_tfestA0 <- tra_tfestA[rowSums(tra_tfestA == 1) >= 3, ]
#tra_tfestA0[tra_tfestA0 > 0] <- 1

filtered_df <- tra_tfestA
#filtered_df[filtered_df > 0] <- 1
filtered_df <- sapply(filtered_df, as.numeric)
sums <- colSums(filtered_df, na.rm = TRUE)
sum_df <- data.frame(columns = colnames(filtered_df), sums = sums)

plot <- ggplot(sum_df, aes(x = columns, y = sums, color = columns)) +
  geom_point(size = 3) +
  labs(title = "Column sums", x = "motifs", y = "sums") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save the plot as a PDF file
ggsave(file=paste0(file_path_out,
                   "sum-mima.pdf"), plot, width = 8, height = 6)


#######################################################################################

sums <- colSums(filtered_df, na.rm = TRUE)
averages <- sums / nrow(filtered_df)

sum_df <- data.frame(columns = colnames(filtered_df), sums = sums, averages = averages)

# Create the barplot using ggplot2
plot2 <- ggplot(sum_df, aes(x = columns, y = averages, color = columns)) +
  geom_point(size = 3) +
  labs(title = "Column Averages", x = "Columns", y = "Average") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save the plot as a PDF file
ggsave(file=paste0(file_path_out,
                   "average-mima.pdf"), plot2, width = 8, height = 6)
#######################################################################################
#head(tra_tfestA0, 10)
#######################################################################################
#######################################################################################
# tra_tfestA0 CONTAINS A RAW SET OF ALL MOTIF MACTHES ACROSS THE DIFFERENT GENES AND CAN BE USED TO FIND BEST COMBINATIONS
#######################################################################################
#######################################################################################

#######################################################################################
#tra_tfestA0$motifs <- apply(tra_tfestA0, 1, function(row) {
#  names(row)[row == 1] %>% paste(collapse = "-")
#})
#fmo_comb <- data.frame(row.names = row.names(tra_tfestA0), motifs = tra_tfestA0$motifs)
#fmo_comb$loc_ID <- row.names(fmo_comb)
#contingency_table_motifs <- table(fmo_comb$motifs, fmo_comb$loc_ID)
#testF =as.data.frame(contingency_table_motifs)
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
testE =as.data.frame(contingency_table_pred_perf)
tfestE = testE %>%
  pivot_wider(names_from = Var2, values_from = Freq)
################################################################################
tfestCDE <- merge( tfestC, 
                   merge(tfestD, tfestE, by = "Var1"), by = "Var1")
tfestCDE$mo_metacluster <- 0
tfestCDE$mo_metacluster[grep("_p1m", tfestCDE$Var1)] <- 1


tfestCDE$prob_rat01 <- ifelse(
  ifelse(tfestCDE$mo_metacluster == 0, 
                              tfestCDE$prob_class0/tfestCDE$prob_class1-1, 
                              tfestCDE$prob_class1/tfestCDE$prob_class0-1) > 0, 1, 0)

tfestCDE$expr_rat01 <- ifelse(
  ifelse(tfestCDE$mo_metacluster == 0, 
                              tfestCDE$expr_class0/tfestCDE$expr_class1-1, 
                              tfestCDE$expr_class1/tfestCDE$expr_class0-1) > 0, 1, 0)

tfestCDE$reliability <- c(tfestCDE$`TRUE`/(tfestCDE$`FALSE`+tfestCDE$`TRUE`))

tfestCDE$chi_expr_class <- apply(tfestCDE[, c("expr_class0", "expr_class1")], 1, function(row) {
  result <- chisq.test(row, p = c(0.5, 0.5))
  result$p.value
})

tfestCDE$chi_prob_class <- apply(tfestCDE[, c("prob_class0", "prob_class1")], 1, function(row) {
  result <- chisq.test(row, p = c(0.5, 0.5))
  result$p.value
})

tfestCDE$chi_TF <- apply(tfestCDE[, c("FALSE", "TRUE")], 1, function(row) {
  result <- chisq.test(row, p = c(0.5, 0.5))
  result$p.value
})

##################################################################################################

write.csv(tfestCDE, file=paste0(file_path_out,
                                "-mima.csv"), row.names=FALSE)


#tfestCDE$perf_row <- rowMeans(tfestCDE[, c(9, 10, 11)])
#tfestCDE$perf_tot<- colMeans(tfestCDE[, c(9, 10, 11)])
