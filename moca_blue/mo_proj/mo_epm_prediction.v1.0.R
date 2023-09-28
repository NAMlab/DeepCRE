# EPM Prediction TEST #
# Use the presence of EPMs as predictor for low and high gene expression 
# EPMs are generalized representations of seqlets with positive or negative 
# contribution scores in deepLift, associated to high and low gene expression prediction
# In general, genes with a positive sum contribution score are predicted to be highly expressed.
# Vice verse, genes with a negative sum contribution score are predicted to be lowly expressed. 
# Consequently, EPMs (associated to postive or negative contrib. scores/ prediction) 
# they can be used equally for interpretation. 
# Here occurences of EPMs per gene are counted. EPMs of p0 are substracted by p1. If the resulting value is positive,
# genes are predicted to be highly expressed and vice verse. 

# Dr. Simon M. Zumkeller 2023-08-25

######## ######### ######### ######### ######## ########
file1 <- "ArthS0Arth-0e3-cwm-W2_gene_none20230825-q1q9.csv"
file2 <- "Arth_S0_predictions.csv" # file with predicted probabilies for expression 0-1
file4 <- "Arabidopsis_thaliana_TPMs-peleke-etal2023.csv" #" file with measurement of models (TPM, quartile classes 0,1,2)
file3 <- "20230825_ArthS0_contrib_scores.csv" # File with average contributions scores for each EPM
dirpath_1 <- "../../ref_seq"
dirpath_2 <- "./out"
dirpath_3 <- "../../mo_nom/out"
file_path_in_file4 <- file.path(dirpath_1, file4)
file_path_in_file3 <- file.path(dirpath_3, file3)
######## ######### ######### ######### ######## ########
#
#
#
#
######## ######### ######### ######### ######## ########
mm0 <- read.table(
  file.path(
    dirpath_2,
    file1),
  header=TRUE,
  sep=",")
#
#
#
mm1 <- mm0[, c("loc_ID",
                      "motif")]
colnames(mm1)[1]<- "loc_ID"
mm1$p0s <- ifelse(sapply(mm1$motif, grepl, pattern = "p0m"), 1, 0)
mm1$p1s <- ifelse(sapply(mm1$motif, grepl, pattern = "p1m"), 1, 0)
mm2 <- mm1 %>%
  group_by(loc_ID) %>%
  summarize(p0_count = sum(p0s), p1_count = sum(p1s))
mm2$epm_pred_val <- mm2$p0_count- mm2$p1_count
mm2$epm_pred <- ifelse(mm2$epm_pred_val >= 0, 1, 0)
mm3 <- mm2[, c("loc_ID",
               "epm_pred")]
value_counts0 <- table(mm3$epm_pred)   # epm_pred value of 0 stands for low expression and 1 for high
########### ############## ############ #####################
#
#
#
#
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
A0<-merge(model0, mm3, by= "loc_ID")
A0$TF_pred <- ifelse(A0$prob == A0$epm_pred, 1, 0) 
##############################################################
#
#
#
model1 <- read.table(
  file_path_in_file4,
  header=TRUE,
  sep=",")
#
##############################################################
##############################################################
ss_model1 <- model1[, c("gene_id", "true_target")]
colnames(ss_model1)[1]<- "loc_ID"
ss_model1 <- subset(ss_model1, true_target !=  2)                                 # !!! CHANGES DATASET SIZE BY HALF !!! # Were not included within the model training.set
A1<-merge(A0, ss_model1, by= "loc_ID")
A1$TF_expr <- ifelse(A1$true_target == A1$epm_pred, 1, 0) 
####
head(A1)
A1$TP_TN <- ifelse(A1$prob==A1$true_target, 1, 0)
A1s <- A1[A1$prob == 1, ]
####
##############################################################
#
#
#
model3 <- read.table(
  file_path_in_file3,
  header=TRUE,
  sep=",")
#
##############################################################
mcon_m0<- merge(mm1,model3, by ="motif")
mcon_m0$p0c<- ifelse(mcon_m0$p0s==1, mcon_m0$contrib_score_aver, 0)
mcon_m0$p1c<- ifelse(mcon_m0$p1s==1, mcon_m0$contrib_score_aver, 0)
mcon_m1<- mcon_m0[, c("loc_ID","p0c","p1c")]
mcon_m2 <- mcon_m1 %>%
  group_by(loc_ID) %>%
  summarize(p0c_count = sum(p0c), p1c_count = sum(p1c))
mcon_m2$epm_contrib_pred_val <- mcon_m2$p0c_count + mcon_m2$p1c_count
mcon_m2$epm_contrib_pred_class <- ifelse(mcon_m2$epm_contrib_pred_val >= 0, 1, 0)
A2 <- merge(A1, mcon_m2, by="loc_ID")                                             # CHANGE A1 to A1s to check expr class
A2$TF_ecpc_expr <- ifelse(A2$epm_contrib_pred_class == A2$true_target, 1, 0)
A2$TF_ecpc_pred <- ifelse(A2$epm_contrib_pred_class == A2$prob, 1, 0)
##############################################################
value_counts1 <- table(A2$TF_ecpc_expr)
value_counts2 <- table(A2$TF_ecpc_pred)
acc <- value_counts1[2]/(value_counts1[2]+value_counts1[1])
acc1 <- value_counts2[2]/(value_counts2[2]+value_counts2[1])
epm_pred_pred <- table(A2$TF_ecpc_expr)
epm_pred_expr <- table(A2$TF_ecpc_pred)
print(epm_pred_pred)
print(acc)
print(epm_pred_expr)
print(acc1)
