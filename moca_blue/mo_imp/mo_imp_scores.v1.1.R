######################
library(rhdf5)
library(tidyr)
library(ggplot2)
library(ggseqlogo)
#setwd("~/ibg-4/Desktop/Rhome/moca_blue/mo_scores")
###################### Setup for "moca_blue" enviroment
NAME0="rdf5_"
SPEC="Soly"
MODEL="S0" # C0 stand for DeepCistrome version 1 (available at 02-may-2023) Standard conditions
DATE= "20230904"
#######################################################
FILE1= "solanum_scores.h5"
FILE2= "solanum_meta_saliency_info.csv"
#######################################################
dirpath_in1 = "../0MOTIFS/saliency_scores/"
dirpath_out = "./out"
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
file_path1 = file.path(dirpath_in1,FILE1)
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
model_parameter <- read.table(
  file.path(
    dirpath_in1,
    FILE2),
  header=TRUE,
  sep=",")
#######################################################
model_parameter$pred_class <- ifelse(model_parameter$pred_prob >= 0.5, 1, 0)
model_parameter$TRUE_class <- ifelse(model_parameter$pred_class == model_parameter$true_target, TRUE, FALSE)
#######################################################
h5file <- H5Fopen(file_path1, "H5F_ACC_RDONLY")
#h5ls(h5file)
saliency_scores <- h5read(h5file,
                            "contrib_scores")
#######################################################
num_arrays <- dim(saliency_scores)[3]
num_columns <- dim(saliency_scores)[2]
result_matrix <- matrix(0, nrow = num_arrays, ncol = num_columns)
for (i in 1:num_arrays) {
  current_array <- saliency_scores[,,i]
  column_sums <- colSums(current_array)
  result_matrix[i,] <- column_sums
}
imp_scores<-as.data.frame(result_matrix)
##################################################################################
model_parameter$sum_imp_score <- rowSums(imp_scores)
model_parameter$sum_imp_score2 <- ifelse(model_parameter$sum_imp_score >= 0, 1, 0)
model_parameter$TRUE_pred_imp <- ifelse(model_parameter$sum_imp_score2 == model_parameter$true_target, TRUE, FALSE)
gene_imp_scores <- as.data.frame(model_parameter$gene_id)
gene_imp_scores <- cbind(gene_imp_scores, imp_scores)
colnames(gene_imp_scores)[1] <- "loc_ID"
##################################################################################
write.table(model_parameter, file = paste0(DATE,NAME0, SPEC, MODEL,"-imp_score_parameter.csv"), sep = ",", col.names = NA, quote = FALSE)
#write.table(gene_imp_scores, file = paste0(DATE,NAME0, SPEC, MODEL,"-imp_score_array.csv"), sep = ",", col.names = NA, quote = FALSE)

#######################################################

selected_gene <- gene_imp_scores[gene_imp_scores[,1] == "Solyc01g111620.3", ]

#selected_gene <- imp_scores[1619, ]
selected_range <- selected_gene[, 750:1500]
# Convert the selected_range matrix to numeric
selected_range <- as.numeric(selected_range)

# Create a data frame for plotting
long_data <- data.frame(
  position = 750:1500,  # Assuming positions start from 1000
  importance = selected_range
)

# Create and display the line plot


line_plot <- ggplot(long_data, aes(x = position, y = importance)) +
  geom_line(size = 0.1, color = "black") +
  geom_vline(xintercept = 1500, color = "red", size = 0.5) +  # Added vertical line
  labs(title = "Range contrib scores",
       x = "Position Index",
       y = "Importance Value") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(750, 1500, by = 50)) +  # Set the x-axis breaks
  theme(panel.grid = element_blank()) +  # Remove grid lines
  scale_y_continuous(limits = c(-0.1, 0.1)) +  # Set y-axis limits
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed")  # Add gray line at y = 0

print(line_plot)

####################

####################
row_numbers <- which(gene_imp_scores[, 1] == "Solyc01g111620.3")

desired_matrix_subset <- saliency_scores[1:4, 750:1500, row_numbers]
head(desired_matrix_subset, 10)
#str(desired_matrix_subset)
rownames(desired_matrix_subset) <- c("A", "C", "G", "T")
####### CONTINUE HERE #########
seq_plot <- ggseqlogo(desired_matrix_subset, method = 'custom', seq_type = 'dna') +
  ylab('1') 


print(seq_plot)
####################
