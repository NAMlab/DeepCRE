# mo_genotype_variance.v1.0.R counts EPM matches per gene across different genotypes.
# This script does not apply any filters like mo_feat-filter.
# Before starting with this tool few preapration had to be made.
# I created a list from differentially and uniformously expressed genes and splitted gene_model.gff files with their annotation: 
# while read id; do grep "$id" *_gene_models.gff >> Sol_genotypes_differentially_expressed_location2.txt; done < Sol_genotypes_differentially_expressed.txt
# For motif search with BLAMM, I created a multifasta file using the extract_range_to_fasta.sh 
# I then merged all fastas into on file and used BLAMM to find epms. 
# Hereafter, I used occ_filter_v1.1.R to parse the matches for further analyses. 
# SMZ 2023-08-30
##############################################################
PROJECT <- "Soly14sols-0e4-cwm-W2"
SPEC <- "Soly"
MODEL <- "S0"
DATE <- "20230830"
##############################################################
word_size <- 14  # This is not a real wordsize like in BLAST but it works similar. All matches smaller than this size will be removed.
# # # # #
file4 <- "outSolS0M0-14sols_0e4-cwm-20230830"    #
file3 <- "Sol_genotypes_differentially_expressed_location2.txt" # 36,6 MB !
file2 <- "Sol_genotypes_uniformally_expressed_location_smaller.txt" # 366 MB !
# # # # # 
dirpath_1 <- "./../../ref_seq/"
dirpath_2 <- "./out"
##############################################################
#   #   #   #   #  #  #   #   #   #  #  #  #  #  #  #  #  #
file_path_in_file3 <- file.path(dirpath_1, file3)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#gene_annot <- import(file_path_in_file3)
gene_annot_diff_expr <- read.table(file_path_in_file3,
                                   header=FALSE,
                                   sep="\t")
split_columns <- strsplit(gene_annot_diff_expr$V1, ":")
split_df <- as.data.frame(do.call(rbind, split_columns))
colnames(split_df) <- paste("Column", 1:ncol(split_df), sep = "")
gene_annot_diff_expr0 <- cbind(split_df, gene_annot_diff_expr[-1])
gene_annot_diff_expr0 <- gene_annot_diff_expr0 %>%
  filter(V3 == "gene")
gene_annot_diff_expr0 <- gene_annot_diff_expr0 %>%
  mutate(gene_id = sub("^[^:]+:([^:]+).*", "\\1", V9))
gene_annot_diff_expr0<- gene_annot_diff_expr0 %>%
  mutate(gene_id = sub("^.*?=(.*?);.*", "\\1", gene_id)) #%>%
#  mutate(gene_id = gsub("ID=", "", gene_id))
gene_annot_diff_expr1 <- gene_annot_diff_expr0[, c(2,5,6,11)]
gene_annot_diff_expr1 <- gene_annot_diff_expr1 %>%
  mutate(V4 = V4 - 1000,
         V5 = V5 + 1000)
gene_annot_diff_expr1$new_column <- paste(gene_annot_diff_expr1$Column2, gene_annot_diff_expr1$V4, sep = ":")
gene_annot_diff_expr1$loc <- paste(gene_annot_diff_expr1$new_column, gene_annot_diff_expr1$V5, sep = "-")
gene_annot_diff_expr1$loc <- gsub(" ", "-", gene_annot_diff_expr1$loc)
gene_annot_diff_expr_locID <- data.frame(loc = gene_annot_diff_expr1$loc, gene_id = gene_annot_diff_expr1$gene_id)
gene_annot_diff_expr_locID <- unique(gene_annot_diff_expr_locID)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##############################################################
#   #   #   #   #  #  #   #   #   #  #  #  #  #  #  #  #  #
file_path_in_file2 <- file.path(dirpath_1, file2)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#gene_annot <- import(file_path_in_file3)
gene_annot_unif_expr <- read.table(file_path_in_file2,
                                   header=FALSE,
                                   sep="\t")
split_columns2 <- strsplit(gene_annot_unif_expr$V1, ":")
split_df2 <- as.data.frame(do.call(rbind, split_columns2))
colnames(split_df2) <- paste("Column", 1:ncol(split_df2), sep = "")
gene_annot_unif_expr0 <- cbind(split_df2, gene_annot_unif_expr[-1])
gene_annot_unif_expr0 <- gene_annot_unif_expr0 %>%
  filter(V3 == "gene")
gene_annot_unif_expr0 <- gene_annot_unif_expr0 %>%
  mutate(gene_id = sub("^[^:]+:([^:]+).*", "\\1", V9))
gene_annot_unif_expr0<- gene_annot_unif_expr0 %>%
  mutate(gene_id = sub("^.*?=(.*?);.*", "\\1", gene_id)) #%>%
#  mutate(gene_id = gsub("ID=", "", gene_id))
gene_annot_unif_expr1 <- gene_annot_unif_expr0[, c(2,5,6,11)]
gene_annot_unif_expr1 <- gene_annot_unif_expr1 %>%
  mutate(V4 = V4 - 1000,
         V5 = V5 + 1000)
gene_annot_unif_expr1$new_column <- paste(gene_annot_unif_expr1$Column2, gene_annot_unif_expr1$V4, sep = ":")
gene_annot_unif_expr1$loc <- paste(gene_annot_unif_expr1$new_column, gene_annot_unif_expr1$V5, sep = "-")
gene_annot_unif_expr1$loc <- gsub(" ", "-", gene_annot_unif_expr1$loc)
gene_annot_unif_expr_locID <- data.frame(loc = gene_annot_unif_expr1$loc, gene_id = gene_annot_unif_expr1$gene_id)
gene_annot_unif_expr_locID <- unique(gene_annot_unif_expr_locID)
# # # # # # # # # # # # # # # # # # # # # # # # # # # 
motif_gene_matches <- read.table(
  file.path(
    dirpath_2,
    file4),
  header=TRUE,
  sep="\t")
motif_gene_matches0 <- motif_gene_matches[motif_gene_matches$mstart <= 1500 |
                                          abs(
                                            motif_gene_matches$mstart-(
                                              motif_gene_matches$gene_end - motif_gene_matches$gene_start)+1) <= 1500, ]
motif_gene_matches0$m_len <- abs(motif_gene_matches0$mstart-motif_gene_matches0$mend)
motif_gene_matches0 <- motif_gene_matches0 %>%
  filter(m_len >= 13)
####                ####                   ####                 ####           ####      ####
merged_diff_df <- merge(motif_gene_matches0,
                   gene_annot_diff_expr_locID,
                   by = "loc")
genomes_diff <- merged_diff_df %>%
  mutate(loc_prefix = sub("^(.*?)_.*", "\\1", loc))
genomes_diff_num <- length(unique(genomes_diff$loc_prefix))
motif_diff_counts <- merged_diff_df %>%
  group_by(gene_id) %>%
  count(motif)
motif_diff_counts0<- motif_diff_counts %>%
  mutate(variance = ifelse(n/genomes_diff_num == floor(n/genomes_diff_num), "conserved", "mutated"))
variance_diff_counts <- table(motif_diff_counts0$variance)
####                ####                   ####                 ####           ####      ####
merged_unif_df <- merge(motif_gene_matches0,
                        gene_annot_unif_expr_locID,
                        by = "loc")
genomes_unif <- merged_unif_df %>%
  mutate(loc_prefix = sub("^(.*?)_.*", "\\1", loc))
genomes_unif_num <- length(unique(genomes_unif$loc_prefix))
motif_unif_counts <- merged_unif_df %>%
  group_by(gene_id) %>%
  count(motif)
motif_unif_counts0<- motif_unif_counts %>%
  mutate(variance = ifelse(n/genomes_unif_num == floor(n/genomes_unif_num), "conserved", "mutated"))
motif_unif_counts0$index <- c(1:nrow(motif_unif_counts0))
motif_unif_counts0
ungrouped_data <- motif_unif_counts0 %>%
  ungroup()
conserved_counts <- numeric(100)
mutated_counts <- numeric(100)
for (i in 1:100) {
  iteration_seed <- sample.int(10^5, 1)
  set.seed(iteration_seed)
  random_subset <- ungrouped_data %>%
    sample_n(size = 1578, replace = FALSE)
  variance_unif_counts <- table(random_subset$variance)
  conserved_counts[i] <- variance_unif_counts["conserved"]
  mutated_counts[i] <- variance_unif_counts["mutated"]
}
average_conserved <- mean(conserved_counts)
average_mutated <- mean(mutated_counts)
summary_table <- c(
  conserved = average_conserved,
  mutated = average_mutated
)
variance_unif_counts <- as.table(summary_table)
# # # # # # # # # # # # # # # # # # # # # # # # # 
variance_diff_df <- as.data.frame(variance_diff_counts)
variance_unif_df <- as.data.frame(variance_unif_counts)
# Merge the two data frames based on the row names
combined_df <- merge(variance_diff_df, variance_unif_df, by = "Var1", all = TRUE)
colnames(combined_df) <- c("Variance", "Diff_Count", "Unif_Count")
combined_df
combined_df$Normalized_Diff_Count <- combined_df$Diff_Count / sum(combined_df$Diff_Count)
combined_df$Normalized_Unif_Count <- combined_df$Unif_Count / sum(combined_df$Unif_Count)
rearranged_df <- combined_df %>%
  pivot_longer(cols = starts_with("Normalized"), 
               names_to = c(".value", "Source"),
               names_sep = "_") %>%
  arrange(Variance, Source)
rearranged_df
### ### ### ### ### ### ### ### ### 
barplot <- ggplot(rearranged_df, aes(x = Source, y = Normalized, fill = Variance)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Source", y = "Percentage", fill = "Variance") +
  ggtitle("Percentages EPM variants amongst uniformously and differentialy expressed genes in Solanum genotypes") +
  scale_fill_manual(values = c("conserved" = "lightgray", "mutated" = "darkgray")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
print(barplot)
####################

head(motif_diff_counts0)
####################
####################

####################

## ## ## ## ## ## ##

####################

####################
####################

# Initialize variables to store results
average_a <- numeric(1000)
average_b <- numeric(1000)
error_rates <- numeric(1000)
# Initialize lists to store data for boxplots
data_for_boxplots_diff <- list()
for (i in 1:1000) {
  # Sample 100 rows randomly
  sampled_indices <- sample(nrow(motif_diff_counts0), 100)
  sampled_motif_diff_counts0s <- motif_diff_counts0[sampled_indices, ]
  
  motif_diff_counts_mutated <- sampled_motif_diff_counts0s[sampled_motif_diff_counts0s$variance == "mutated", ]
  motif_diff_counts_conserved <- sampled_motif_diff_counts0s[sampled_motif_diff_counts0s$variance == "conserved", ]
  genes_conserved_with_diff_expr <- anti_join(motif_diff_counts_conserved, motif_diff_counts_mutated, by = "gene_id")
  genes_mutated_with_diff_expr <- anti_join(motif_diff_counts_mutated, genes_conserved_with_diff_expr, by = "gene_id")
  a <- length(unique(genes_conserved_with_diff_expr$gene_id))
  b <- length(unique(genes_mutated_with_diff_expr$gene_id))
  average_a[i] <- a
  average_b[i] <- b
  error_rates[i] <- b / (a + b)
  
  # Store data for boxplot
  data_for_boxplots_diff[[i]] <- c(a, b)
}

# Calculate average and error rate
avg_a <- mean(average_a)
avg_b <- mean(average_b)
avg_error_rate_diff <- mean(error_rates)
error_rate_sd_diff <- sd(error_rates)
# Calculate error rates for averages of a and b
error_rate_avg_a <- avg_b / (avg_a + avg_b)
error_rate_avg_b <- avg_a / (avg_a + avg_b)
# Print results
cat("Average a:", avg_a, "\n")
cat("Average b:", avg_b, "\n")
cat("Average Error Rate Diff:", avg_error_rate_diff, "\n")
cat("Error Rate Standard Deviation:", error_rate_sd_diff, "\n")
# Print error rates for averages of c and d
cat("Error Rate for Average of c:", error_rate_avg_a, "\n")
cat("Error Rate for Average of d:", error_rate_avg_b, "\n")

# Create a dataframe from your data
results_df0 <- data.frame(
  Category = c("Average_a", "Average_b"),
  Average = c(avg_a, avg_b),
  Standard_Deviation = c(error_rate_sd_diff, error_rate_sd_diff)
)
results_df0

head(data_for_boxplots_diff)

# Calculate statistics for 'a' (left values) and 'b' (right values)
a_values <- sapply(data_for_boxplots_diff, "[[", 1)
b_values <- sapply(data_for_boxplots_diff, "[[", 2)

# Calculate statistics for 'a'
a_min <- min(a_values)
a_max <- max(a_values)
a_mean <- mean(a_values)
a_q1 <- quantile(a_values, 0.25)
a_q3 <- quantile(a_values, 0.75)

# Calculate statistics for 'b'
b_min <- min(b_values)
b_max <- max(b_values)
b_mean <- mean(b_values)
b_q1 <- quantile(b_values, 0.25)
b_q3 <- quantile(b_values, 0.75)

# Print the results
cat("Statistics for 'a' (left values):\n")
cat("Min:", a_min, "\n")
cat("Max:", a_max, "\n")
cat("Mean:", a_mean, "\n")
cat("Q1 (25th percentile):", a_q1, "\n")
cat("Q3 (75th percentile):", a_q3, "\n")

cat("\nStatistics for 'b' (right values):\n")
cat("Min:", b_min, "\n")
cat("Max:", b_max, "\n")
cat("Mean:", b_mean, "\n")
cat("Q1 (25th percentile):", b_q1, "\n")
cat("Q3 (75th percentile):", b_q3, "\n")
####

####

####
####################

# Initialize variables to store results
average_b <- numeric(1000)
average_d <- numeric(1000)
error_rates <- numeric(1000)
# Initialize lists to store data for boxplots
data_for_boxplots_unif <- list()
# Perform the process 1000 times
for (i in 1:1000) {
  # Sample 100 rows randomly
  sampled_indices <- sample(nrow(motif_unif_counts0), 100)
  sampled_motif_unif_counts0s <- motif_unif_counts0[sampled_indices, ]
  
  motif_unif_counts_mutated <- sampled_motif_unif_counts0s[sampled_motif_unif_counts0s$variance == "mutated", ]
  motif_unif_counts_conserved <- sampled_motif_unif_counts0s[sampled_motif_unif_counts0s$variance == "conserved", ]
  genes_conserved_with_unif_expr <- anti_join(motif_unif_counts_conserved, motif_unif_counts_mutated, by = "gene_id")
  genes_mutated_with_unif_expr <- anti_join(motif_unif_counts_mutated, genes_conserved_with_unif_expr, by = "gene_id")
  a <- length(unique(genes_conserved_with_unif_expr$gene_id))
  b <- length(unique(genes_mutated_with_unif_expr$gene_id))
  average_a[i] <- a
  average_b[i] <- b
  error_rates[i] <- b / (a + b)
  
  # Store data for boxplot
  data_for_boxplots_unif[[i]] <- c(a, b)
}


# Calculate average and error rate
avg_a <- mean(average_a)
avg_b <- mean(average_b)
avg_error_rate_unif <- mean(error_rates)
error_rate_sd_unif <- sd(error_rates)
# Calculate error rates for averages of a and b
error_rate_avg_a <- avg_b / (avg_a + avg_b)
error_rate_avg_b <- avg_a / (avg_a + avg_b)
# Print results
cat("Average a:", avg_a, "\n")
cat("Average b:", avg_b, "\n")
cat("Average Error Rate unif:", avg_error_rate_unif, "\n")
cat("Error Rate Standard Deviation:", error_rate_sd_unif, "\n")
# Print error rates for averages of c and d
cat("Error Rate for Average of c:", error_rate_avg_a, "\n")
cat("Error Rate for Average of d:", error_rate_avg_b, "\n")

# Create a dataframe from your data
results_df0 <- data.frame(
  Category = c("Average_a", "Average_b"),
  Average = c(avg_a, avg_b),
  Standard_Deviation = c(error_rate_sd_unif, error_rate_sd_unif)
)
results_df0

head(data_for_boxplots_unif)

# Calculate statistics for 'a' (left values) and 'b' (right values)
c_values <- sapply(data_for_boxplots_unif, "[[", 1)
d_values <- sapply(data_for_boxplots_unif, "[[", 2)

# Calculate statistics for 'a'
c_min <- min(c_values)
c_max <- max(c_values)
c_mean <- mean(c_values)
c_q1 <- quantile(c_values, 0.25)
c_q3 <- quantile(c_values, 0.75)

# Calculate statistics for 'b'
d_min <- min(d_values)
d_max <- max(d_values)
d_mean <- mean(d_values)
d_q1 <- quantile(d_values, 0.25)
d_q3 <- quantile(d_values, 0.75)

# Print the results
cat("Statistics for 'a' (left values):\n")
cat("Min:", c_min, "\n")
cat("Max:", c_max, "\n")
cat("Mean:", c_mean, "\n")
cat("Q1 (25th percentile):", c_q1, "\n")
cat("Q3 (75th percentile):", c_q3, "\n")

cat("\nStatistics for 'b' (right values):\n")
cat("Min:", d_min, "\n")
cat("Max:", d_max, "\n")
cat("Mean:", d_mean, "\n")
cat("Q1 (25th percentile):", d_q1, "\n")
cat("Q3 (75th percentile):", d_q3, "\n")

# Create a data frame for boxplot
boxplot_data_diff <- data.frame(
  Category = factor(c("genes_diff_expr_EPM_conserved", "genes_diff_expr_EPM_mutated")),
  Min = c(a_min, b_min),
  Q1 = c(a_q1, b_q1),
  Median = c(a_mean, b_mean),
  Q3 = c(a_q3, b_q3),
  Max = c(a_max, b_max)
)
#####

# Create a data frame for boxplot
boxplot_data_unif <- data.frame(
  Category = factor(c("genes_unif_expr_EPM_conserved", "genes_unif_expr_EPM_mutated")),
  Min = c(c_min, d_min),
  Q1 = c(c_q1, d_q1),
  Median = c(c_mean, d_mean),
  Q3 = c(c_q3, d_q3),
  Max = c(c_max, d_max)
)
# Combine the two data frames into one
combined_boxplot_data <- bind_rows(
  boxplot_data_unif,
  boxplot_data_diff
)
# Define custom colors for 'a', 'b', 'c', and 'd'
custom_colors <- c("genes_diff_expr_EPM_conserved" = "cyan3",
                   "genes_diff_expr_EPM_mutated" = "gold",
                   "genes_unif_expr_EPM_conserved" = "cyan4",
                   "genes_unif_expr_EPM_mutated" = "darkgoldenrod2")

ggplot(combined_boxplot_data, aes(x = Category, ymin = Min, lower = Q1, middle = Median, upper = Q3, ymax = Max, fill = Category)) +
  geom_boxplot(stat = "identity", width = 0.5, position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = custom_colors) +  # Assign custom colors
  labs(x = "Category", y = "Value") +
  ggtitle("Occurences of EPM states for uniform and different gene expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major.x = element_blank())  # Remove vertical gridlines

################################################################################################

motif_diff_counts_mutated0 <- motif_diff_counts0[motif_diff_counts0$variance == "mutated", ]
motif_diff_counts_conserved0 <- motif_diff_counts0[motif_diff_counts0$variance == "conserved", ]
genes_conserved_with_diff_expr0 <- anti_join(motif_diff_counts_conserved0, motif_diff_counts_mutated0, by = "gene_id")
genes_mutated_with_diff_expr0 <- anti_join(motif_diff_counts_mutated0, genes_conserved_with_diff_expr0, by = "gene_id")

motif_unif_counts_mutated0 <- motif_unif_counts0[motif_unif_counts0$variance == "mutated", ]
motif_unif_counts_conserved0 <- motif_unif_counts0[motif_unif_counts0$variance == "conserved", ]
genes_conserved_with_unif_expr0 <- anti_join(motif_unif_counts_conserved0, motif_unif_counts_mutated0, by = "gene_id")
genes_mutated_with_unif_expr0 <- anti_join(motif_unif_counts_mutated0, genes_conserved_with_unif_expr0, by = "gene_id")

# Add a column indicating the source data frame
genes_conserved_with_diff_expr0$source_df <- "genes_conserved_with_diff_expr0"
genes_mutated_with_diff_expr0$source_df <- "genes_mutated_with_diff_expr0"
genes_conserved_with_unif_expr0$source_df <- "genes_conserved_with_unif_expr0"
genes_mutated_with_unif_expr0$source_df <- "genes_mutated_with_unif_expr0"

# Combine the data frames into one
combined_genes_df <- bind_rows(
  genes_conserved_with_diff_expr0,
  genes_mutated_with_diff_expr0,
  genes_conserved_with_unif_expr0,
  genes_mutated_with_unif_expr0
)
head(combined_genes_df)
write.csv(combined_genes_df, "combined_genes_data.csv", row.names = FALSE)
