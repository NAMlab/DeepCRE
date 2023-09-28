#####################
# Extract the contribution scores of EPMs
# from hdf5 files from modisco
# Dr. SM Zumkeller 2023-08-25
######################


#setwd("~/ibg-4/Desktop/Rhome/moca_blue/mo_nom")
###################### Setup for "moca_blue" enviroment
NAME0="rdf5_CWMs"
SPEC="Arth"
MODEL="M0" # C0 stand for DeepCistrome version 1 (available at 02-may-2023) Standard conditions
DATE = "20230904"
#
#
#
library(dplyr)
library(rhdf5)
#######################################################
FILE1= "Arabidopsis_MSR_modisco.hdf5"
#######################################################
dirpath_in = "../0MOTIFS/MODISCO_MSR"
dirpath_out = "./out"
#
file_path_out <- file.path(dirpath_out, paste0(DATE,"_",SPEC,MODEL,"_contrib_scores"))
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
h5file <- H5Fopen(file.path(
  dirpath_in,
  FILE1), "H5F_ACC_RDONLY")
#h5ls(h5file)
metacluster_group <- h5read(h5file,
                            "metacluster_idx_to_submetacluster_results")
#######################################################
# loop through the metaclusters 0 and 1
for (i in names(metacluster_group)) {
  metacluster <- metacluster_group[[i]]
  patterns = metacluster[['seqlets_to_patterns_result']][['patterns']]
}
# Define the pattern names to iterate over
########################################################
length(patterns[['all_pattern_names']])
x0 = length(metacluster_group[["metacluster_0"]][["seqlets_to_patterns_result"]][["patterns"]])-1
x1 = length(metacluster_group[["metacluster_1"]][["seqlets_to_patterns_result"]][["patterns"]])-1
y01 = x0+x1
y02 = x0*2+x1*2
###################################################
pattern_names <- paste0("pattern_", 0:y01) ############################################## !!! MANUAL ADJ REQUIRED
# Initialize an empty list to store the matrices
matricesF0 <- list()
matricesF1 <- list()
matricesR0 <- list()
matricesR1 <- list()
# Loop over the pattern names
for (pattern_name in pattern_names) {
  # Extract the matrix for the current pattern name
  matrixF0 <- metacluster_group[["metacluster_0"]][["seqlets_to_patterns_result"]][["patterns"]][[pattern_name]][["task0_contrib_scores"]][["fwd"]]
  matrixF1 <- metacluster_group[["metacluster_1"]][["seqlets_to_patterns_result"]][["patterns"]][[pattern_name]][["task0_contrib_scores"]][["fwd"]]
  matrixR0 <- metacluster_group[["metacluster_0"]][["seqlets_to_patterns_result"]][["patterns"]][[pattern_name]][["task0_contrib_scores"]][["rev"]]
  matrixR1 <- metacluster_group[["metacluster_1"]][["seqlets_to_patterns_result"]][["patterns"]][[pattern_name]][["task0_contrib_scores"]][["rev"]]
  #pfms in sequence
  #cwms in sequence [task0_contrib_scores] : contibution weight matrices
  # Append the matrix to the list
  matricesF0[[pattern_name]] <- matrixF0
  matricesF1[[pattern_name]] <- matrixF1
  matricesR0[[pattern_name]] <- matrixR0
  matricesR1[[pattern_name]] <- matrixR1
}
##################################################################################### 
###############                ASSIGN NOMENCLATURE            ###################
###############                ASSIGN NOMENCLATURE            ###################
###############                ASSIGN NOMENCLATURE            ###################

##################################################################################### PFM to PWM to CWM
seqletls_lengths_p1 <- list()
# Loop through each pattern name and get the length of seqletls
for (pattern_name in pattern_names) {
  seqletls <- metacluster_group[["metacluster_1"]][["seqlets_to_patterns_result"]][["patterns"]][[pattern_name]][["seqlets_and_alnmts"]][["seqlets"]]
  seqletls_lengths_p1[[pattern_name]] <- length(seqletls)
}

seqletls_lengths_p0 <- list()
# Loop through each pattern name and get the length of seqletls
for (pattern_name in pattern_names) {
  seqletls <- metacluster_group[["metacluster_0"]][["seqlets_to_patterns_result"]][["patterns"]][[pattern_name]][["seqlets_and_alnmts"]][["seqlets"]]
  seqletls_lengths_p0[[pattern_name]] <- length(seqletls)
}

################
names(matricesF0) <- paste0(names(matricesF0),
                            "_p0m",
                            sprintf("%02d",
                                    as.numeric(substring(names(matricesF0), 9))))
names(matricesF0) <- substr(names(matricesF0), nchar(names(matricesF0)) - 4, nchar(names(matricesF0)))
names(matricesF0) <- paste0(names(matricesF0),
                            "F")
names(matricesF0) <- paste0("epm_", SPEC, "_", MODEL, "_", names(matricesF0))
###############                  ####################            ###################
names(matricesF1) <- paste0(names(matricesF1),
                            "_p1m",
                            sprintf("%02d",
                                    as.numeric(substring(names(matricesF1), 9))))
names(matricesF1) <- substr(names(matricesF1), nchar(names(matricesF1)) - 4, nchar(names(matricesF1)))
names(matricesF1) <- paste0(names(matricesF1),
                            "F")
names(matricesF1) <- paste0("epm_", SPEC, "_", MODEL, "_", names(matricesF1))
###############                  ####################            ###################
names(matricesR1) <- paste0(names(matricesR1),
                            "_p1m",
                            sprintf("%02d",
                                    as.numeric(substring(names(matricesR1), 9))))
names(matricesR1) <- substr(names(matricesR1), nchar(names(matricesR1)) - 4, nchar(names(matricesR1)))
names(matricesR1) <- paste0(names(matricesR1),
                            "R")
names(matricesR1) <- paste0("epm_", SPEC, "_", MODEL, "_", names(matricesR1))
###############                  ####################            ###################
names(matricesR0) <- paste0(names(matricesR0),
                            "_p0m",
                            sprintf("%02d",
                                    as.numeric(substring(names(matricesR0), 9))))
names(matricesR0) <- substr(names(matricesR0), nchar(names(matricesR0)) - 4, nchar(names(matricesR0)))
names(matricesR0) <- paste0(names(matricesR0),
                            "R")
names(matricesR0) <- paste0("epm_", SPEC, "_", MODEL, "_", names(matricesR0))
###############                  ####################            ###################   # ADD NUMBER OF SEQLETS TO NAMES !!!!
seqlets_count_p0 <- head(seqletls_lengths_p0, x0)

for (i in seq_along(seqlets_count_p0)) {
  name <- paste0(names(matricesF0)[i], "_", seqlets_count_p0[[i]])
  names(matricesF0)[i] <- name
}

for (i in seq_along(seqlets_count_p0)) {
  name <- paste0(names(matricesR0)[i], "_", seqlets_count_p0[[i]])
  names(matricesR0)[i] <- name
}

seqlets_count_p1 <- head(seqletls_lengths_p1, x1)

for (i in seq_along(seqlets_count_p1)) {
  name <- paste0(names(matricesF1)[i], "_", seqlets_count_p1[[i]])
  names(matricesF1)[i] <- name
}

for (i in seq_along(seqlets_count_p1)) {
  name <- paste0(names(matricesR1)[i], "_", seqlets_count_p1[[i]])
  names(matricesR1)[i] <- name
}
#################################################################################### contribution scores, to weitgh matrix
motifs <- c(matricesF0,matricesF1,matricesR0,matricesR1)
#####################################################################################
results_df1 <- data.frame(matrix_name = character(0), total_sum = numeric(0))
for (matrix_name in names(matricesF0)) {
  matrix <- matricesF0[[matrix_name]]
  row_sums <- rowSums(matrix)
  col_sums <- colSums(matrix)
  total_sum <- 1*(sum(row_sums) + sum(col_sums))  # Calculate the total sum
  max_val <- max(matrix)
  min_val <- min(matrix)
  # Append the result to the data frame
  results_df1 <- rbind(results_df1, data.frame(motif = matrix_name,
                                             contrib_score_sum = total_sum,
                                             contrib_score_max = max_val,
                                             contrib_score_min = min_val))
}
contrib_scores_F0 <- results_df1
#####################################################################################
results_df2 <- data.frame(matrix_name = character(0), total_sum = numeric(0))
for (matrix_name in names(matricesF1)) {
  matrix <- matricesF1[[matrix_name]]
  row_sums <- rowSums(matrix)
  col_sums <- colSums(matrix)
  total_sum <- 1*((sum(row_sums) + sum(col_sums)))  # Calculate the total sum
  max_val <- max(matrix)
  min_val <- min(matrix)
  # Append the result to the data frame
  results_df2 <- rbind(results_df2, data.frame(motif = matrix_name,
                                             contrib_score_sum = total_sum,
                                             contrib_score_max = max_val,
                                             contrib_score_min = min_val))
}
contrib_scores_F1 <- results_df2
#####################################################################################
results_df3 <- data.frame(matrix_name = character(0), total_sum = numeric(0))
for (matrix_name in names(matricesR0)) {
  matrix <- matricesR0[[matrix_name]]
  row_sums <- rowSums(matrix)
  col_sums <- colSums(matrix)
  total_sum <- ((sum(row_sums) + sum(col_sums)))  # Calculate the total sum
  max_val <- max(matrix)
  min_val <- min(matrix)
  # Append the result to the data frame
  results_df3 <- rbind(results_df3, data.frame(motif = matrix_name,
                                             contrib_score_sum = total_sum,
                                             contrib_score_max = max_val,
                                             contrib_score_min = min_val))
}
contrib_scores_R0 <- results_df3
#####################################################################################
results_df4 <- data.frame(matrix_name = character(0), total_sum = numeric(0))
for (matrix_name in names(matricesR1)) {
  matrix <- matricesR1[[matrix_name]]
  row_sums <- rowSums(matrix)
  col_sums <- colSums(matrix)
  total_sum <- 1*((sum(row_sums) + sum(col_sums)))  # Calculate the total sum
  max_val <- max(matrix)
  min_val <- min(matrix)
  # Append the result to the data frame
  results_df4 <- rbind(results_df4, data.frame(motif = matrix_name,
                                             contrib_score_sum = total_sum,
                                             contrib_score_max = max_val,
                                             contrib_score_min = min_val))
}
contrib_scores_R1 <- results_df4
#####################################################################################
contrib_score_table <- rbind(contrib_scores_F0, contrib_scores_R0, contrib_scores_F1, contrib_scores_R1)



write.csv(contrib_score_table, file=paste0(file_path_out,
                                ".csv"), row.names=FALSE)