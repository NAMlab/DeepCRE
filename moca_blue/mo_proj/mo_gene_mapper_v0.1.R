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
file1 <- "SolyMSR_on_Spe-ch01-0e3_gene_none20230530-q1q9.csv"
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
library(reshape2)


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

################################################################################
a0_fd <- merg_df01mm[, c(5,1,2,3,4,10,11,6,12)]

# Determine the percentile cutoff values
lower_percentile <- quantile(a0_fd$prob, 0.25)
upper_percentile <- quantile(a0_fd$prob, 0.75)

# Subset rows based on percentiles
a0_fd_ssprob <- a0_fd[a0_fd$prob<= lower_percentile | a0_fd$prob >= upper_percentile, ]
# Determine the percentile cutoff values
lower_percentile1 <- quantile(a0_fd$logMaxTPM, 0.25)
upper_percentile1 <- quantile(a0_fd$logMaxTPM, 0.75)
# Subset rows based on percentiles
a0_fd_ssTPM <- a0_fd[a0_fd$logMaxTPM<= lower_percentile1 | a0_fd$logMaxTPM >= upper_percentile1, ]

a0_fd_ssprobTPM2 <- a0_fd_ssprob[a0_fd_ssprob$logMaxTPM<= lower_percentile1 |
                                   a0_fd_ssprob$logMaxTPM >= upper_percentile1, ]


###############################################################################
head(a0_fd_ssprobTPM2)
# Extracting values from the "motif" column and shortening them
a0_fd_ssprobTPM2 <- a0_fd_ssprobTPM2 %>% 
  mutate(motif0 = substr(motif, 4, 19))
# Adding "strand.x" and "region" values to "motif0"
a0_fd_ssprobTPM2 <- a0_fd_ssprobTPM2 %>% 
  mutate(motif0 = paste(motif0, strand.x, ifelse(region == "upstream", "up", "do"), sep = ""))
# Creating the "motif1" column by combining "motif0" and "dist_transc_border"
a0_fd_ssprobTPM2 <- a0_fd_ssprobTPM2 %>% 
  mutate(motif1 = paste(motif0, dist_transc_border, sep = "/"))
a1_fd <- a0_fd_ssprobTPM2[, c(2,10,4)]
a1_fd_unique <- distinct(a1_fd)
a1_fd_reshaped <- a1_fd_unique %>%
  pivot_wider(names_from = motif0, values_from = dist_transc_border)  #### THIS IS WORKING SOMEHOW... FRAGILE! HANDLE WITH CARE (GREAT!)
################################################################################

a2<-as.data.frame(a1_fd_reshaped)
rownames(a2)<-(a2$loc_ID)
a2_0 <- as.data.frame(t(a2))
# Remove the first row
a2_0 <- a2_0[-1, ]
# Rename the first column to "motif0"
a2_0$motif0 <- rownames(a2_0)
# Separate "motif0" into three new columns
a2_0 <- cbind(
  a2_0,
  motif = substr(a2_0$motif0, 1, 16),
  strand = ifelse(grepl("\\+", a2_0$motif0), "+", "-"),
  region = ifelse(grepl("up", a2_0$motif0), "up", "do")
)
# Reorder the columns
a2_0 <- a2_0[, c("motif", "strand", "region", colnames(a2_0)[-1])]
# Export the updated DataFrame as CSV
#write.csv(a2_0, file = "mo_dist_matrix_preliminary.csv", row.names = FALSE)
# Remove the row and column used for row and column names
# Coerce columns to character type
a2_fixed <- as.data.frame(lapply(a2_0, as.character), stringsAsFactors = FALSE)

# Replace "NULL" values with zeros
a2_fixed[a2_fixed == "NULL"] <- NA

# Export the fixed DataFrame as CSV
write.csv(a2_fixed, file = "mo_dist_matrix_preliminary.csv", row.names = FALSE)

################################################################################

############################################################################### NOTE! DATASET HAS BEEN FILTERED FOR TRUE PERCENTILES (0.05/0.095!!!!

# Contingency table of "motif" and "loc_ID"
contingency_table_loc_ID <- table( a0_fd_ssprobTPM2$loc_ID, a0_fd_ssprobTPM2$motif)
# Contingency table of "motif" and "GO"
#contingency_table_GO <- table(a0_fd_ssprobTPM2$motif, a0_fd_ssprobTPM2$NAME)
# Contingency table of "motif" and "class"

testA =as.data.frame(contingency_table_loc_ID)
tfestA = testA %>%
  pivot_wider(names_from = Var2, values_from = Freq)

#tfestA[tfestA != 0] <- tfestA$Var1

# Create a copy of the original dataframe
tfestA0 <- tfestA

# Convert all columns except Var1 to character type
tfestA0[, -1] <- lapply(tfestA0[, -1], as.character)

# Iterate over columns and replace non-zero values
for (i in 2:ncol(tfestA0)) {
  tfestA0[tfestA0[, i] != "0", i] <- tfestA0$Var1[tfestA0[, i] != "0"]
}

# View the modified dataframe
tfestA0

write.csv(tfestA, file = "mo_binary_matrix_preliminary.csv", row.names = FALSE)

write.csv(tfestA0, file = "mo_gene_list_preliminary.csv", row.names = FALSE)



###############################################################################################

###############################################################################################
