#####################
# This script extract motifs in CWM format from hdf5 files and stores them in jaspar format.
# The script also names the motifs according to proposed nomenclature and the selected input specs
######################
library(rhdf5)
library(tidyr)
###################### Setup for "moca_blue" enviroment
NAME0="rdf5_epm"
SPEC="Zema"
MODEL="C0" # C0 stand for DeepCistrome version 1 (available at 02-may-2023) Standard conditions
#######################################################
FILE1= "modisco.hdf5"
#######################################################
dirpath_in = "../Motifs/MOTIFS_DC1_ZEA"
dirpath_out = "./out"
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
#paste0(data$motif, "m", sprintf("%02d", as.numeric(substring(data$pattern, 9))))
###############                  ####################            ###################
# Assuming seqletls_lengths_p0 is a list of pattern lengths and matricesF0 is a list of matrices with matching pattern names
#for (pattern_name in names(matricesF0)) {
#  if (pattern_name %in% names(seqletls_lengths_p0)) {
#    matricesF0[[pattern_name]] <- matricesF0[[pattern_name]] * seqletls_lengths_p0[[pattern_name]]
#  }
#}
# Assuming seqletls_lengths_p0 is a list of pattern lengths and matricesF0 is a list of matrices with matching pattern names
#for (pattern_name in names(matricesF1)) {
#  if (pattern_name %in% names(seqletls_lengths_p1)) {
#    matricesF1[[pattern_name]] <- matricesF1[[pattern_name]] * seqletls_lengths_p1[[pattern_name]]
#  }
#}
##############
#for (pattern_name in names(matricesR0)) {
#  if (pattern_name %in% names(seqletls_lengths_p0)) {
#    matricesR0[[pattern_name]] <- matricesR0[[pattern_name]] * seqletls_lengths_p0[[pattern_name]]
#  }
#}
# Assuming seqletls_lengths_p0 is a list of pattern lengths and matricesF0 is a list of matrices with matching pattern names
#for (pattern_name in names(matricesR1)) {
#  if (pattern_name %in% names(seqletls_lengths_p1)) {
#    matricesR1[[pattern_name]] <- matricesR1[[pattern_name]] * seqletls_lengths_p1[[pattern_name]]
#  }
#}
  ###############                ASSIGN NOMENCLATURE            ###################
 ###############                ASSIGN NOMENCLATURE            ###################
###############                ASSIGN NOMENCLATURE            ###################
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

#m0 <- motifs[1]
#m1 <- motifs[[1]]
#name <- names(m0)[1]
#seq_count <- sub(".*_([0-9]+)$", "\\1", name)
#nfcwm <- abs(m1)
#nfcwm <- as.numeric(seq_count)*(nfcwm/max(nfcwm))
####################################################################################
for (i in seq_along(motifs)) {
  m0 <- motifs[i]
  m1 <- motifs[[i]]
  name <- names(m0)[1]
  seq_count <- sub(".*_([0-9]+)$", "\\1", name)
  nfcwm <- abs(m1)
  nfcwm <- round(as.numeric(seq_count)*(nfcwm/max(nfcwm)))
  motifs[[i]] <- abs(nfcwm)
}

####################################################################################
pfms<- array(unlist(motifs),dim = c(4, 14, y02))    # make correction here - not pfm but cwm!!!!! extracting from sequence gives the PPM ! position probability matrix!!!!
ls_pfms<- list()
for (idx in seq(1:y02)){
  ls_pfms[[idx]] <-pfms[, , idx]
}
ls_pfms_str <- lapply(ls_pfms, function(x) {
  apply(x, c(1, 2), as.character)
})
#####################################################################################
create_text <- function(m){
  res <- paste0(">motif", "\n")
  rows <- c('A', 'C', 'G', 'T')
  for(i in 1:nrow(m)){
    res <- paste0(res, paste0(rows[i],' ', paste0('[', paste(as.character(m[i, 1:ncol(m)]), collapse = "\t"), ']', "\n")))
  }
  return(res)
}
##########################                            ###############################
                          ############################
text <- lapply(ls_pfms, create_text)
for (idx in seq(1:y02)){
  text[[idx]] <- gsub("motif", paste0(names(motifs)[idx]), text[[idx]])
}
writeLines(unlist(text), paste0(NAME0,SPEC,MODEL,"_cwm-motifs.jaspar"))
#####################################################################################
