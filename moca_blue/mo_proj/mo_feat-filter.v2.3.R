# This scirpt should run from the directory:
# /home/ibg-4/Desktop/Rhome/solanum_motifs
# There are four files in total which will'be used to process the data
# There is SolyMSR-TSS_motif_position.csv and SolyMSR-TTS_motif_position.csv
# These two files contain the summary statistics of the motifs and their margins
# Then there is Solanum_lycopersicum.SL3.0.55.chr.gff3 (AND do it for spenn_v2.0_gene_models_annot.gff)
# These files contain information about the genes region freature
# These files are located in the directory spenn_v2.0_gene_models_annot.gff
# Finally, there is the file out_all.txt
# this file is located in the directory /home/ibg-4/Desktop/Rhome/solanum_motifs/out
# It contains which motifs matched what gene region
# All files need to imported into R and combined
# The goal of this script is to characterize genes that have identified motifs in their flanking regions

#setwd("/home/ibg-4/Desktop/Rhome/solanum_motifs")
#gff annotations might be edited eg. sed -i 's/ID=gene:[^A]*AT/ID=gene:AT/g' path/to/file
##############################################################
PROJECT <- "ArthS0_Arly_0e3-cwm"
SPEC <- "Arly"
MODEL <- "S0"
DATA <- "20230707"
##############################################################
dirpath_1 <- "../../ref_seq"
dirpath_2 <- "./out"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   #   #   #   #  #  #   #   #   #  #  #  #  #  #  #  #  #
file4 <- "outArthS0-Arly_0e3-cwm.txt"    #   These are the filtered results of the the BLAMM output occurence.txt
#   #   #   #   #  #  #   #   #   #  #  #  #  #  #  #  #  #
file1 <- "ArthS0-TSS_motif_ranges.csv"         # These are the summary statistics of the seqlet distribution of the HDF5 file  
file2 <- "ArthS0-TTS_motif_ranges.csv"         # These are the summary statistics of the seqlet distribution of the HDF5 file
file3 <- "Arabidopsis_lyrata.v.1.0.56.chr-mod.gff3" # Genome annotation file
#   #   #   #   #  #  #   #   #   #  #  #  #  #  #  #  #  #
#   #   #   #   #  #  #   #   #   #  #  #  #  #  #  #  #  #
# Available filters for motif_orient: "forward", "reverse", "none"
Filter_motif_orient <- "none"
# Available filters for annot_type : "gene", "CDS", mRNA", "UTR", "none"
Filter_annot_type   <- "gene"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#ratio_treshold <- 10
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##############################################################
##############################################################
# Load the necessary packages
library(dplyr)
library(stringr)
##############################################################
#THIS MUST BE CHANGED AN ADJUSTED WHEN IT IS USED ON THE OTHER SPECIES !!
convert_epm <- function(code) {
  code1 <- substr(code, 1, 12) 
  code2 <- substr(code, 13, 17)
  code <- paste0(code1, code2)
  
  code <- gsub("Sola", "Soly", code, fixed = TRUE) #ONLY SOLA mispelling
  code <- gsub("epm", "epm_", code, fixed = TRUE)
  code <- gsub("__", "_", code, fixed = TRUE)
  return(code)
}
##############################################################

# Read in the motif position files
tss_motifs <- read.table(file1, header = TRUE
                         , sep=",")
tts_motifs <- read.table(file2, header = TRUE
                         , sep=",")
#tss_motifs <- subset(tss_motifs, ratio > ratio_treshold/100)
#tts_motifs <- subset(tts_motifs, ratio > ratio_treshold/100)
#ratio_treshold
# Read in the gene annotation file
gene_annot <- read.table(
  file.path(
    dirpath_1, file3),
  header=FALSE, sep="\t")
colnames(gene_annot) <- c("chr",
                          "source",
                          "type",
                          "start",
                          "end",
                          "score",
                          "strand",  
                          "phase",
                          "attributes")
gene_annotX <- gene_annot
############################################ This is superflous and can made shorter ##########
if (Filter_annot_type == "gene") {
  gene_annot <- gene_annot %>%
    filter(grepl("gene", type))
} else if (Filter_annot_type == "CDS") {
  gene_annot <- gene_annot %>%
    filter(grepl("CDS", type))
} else if (Filter_annot_type == "mRNA") {
  gene_annot <- gene_annot %>%
    filter(grepl("mRNA", type))
} else if (Filter_annot_type == "UTR") {
  gene_annot <- gene_annot %>%
    filter(grepl("UTR", type))
} else if (Filter_annot_type == "none") {
  gene_annot <- gene_annot %>%
    filter(grepl("gene","CDS","mRNA","UTR", type))
  # No filtering applied, retain all rows
  # You can remove this condition if no filter is needed for "none"
} else {
  # Handle cases where Filter_annot_type has unexpected values
  stop("Invalid value for Filter_annot_type")
}
# Read in the motif-gene matches file
motif_gene_matches <- read.table(
  file.path(
    dirpath_2,
    file4),
  header=TRUE,
  sep="\t")
##############################################################
#Combine DFs by Chr:start-end as loc; Use extracted ranges
#Create loc to idx df
gene_annot$loc <- paste(
  gene_annot$chr,
  ":",
  gene_annot$start - 1000,
  "-",
  gene_annot$end + 1000,
  sep = "")
#Reduce attributes
gene_annot$loc_ID <- str_extract(
  gene_annot$attributes,
  "(?<=ID=).*?(?=;)")
gene_annot$loc_ID <- ifelse(
  grepl(":",
        gene_annot$loc_ID),
  sub(".*:", "",
      gene_annot$loc_ID),
  gene_annot$loc_ID)
#subset the df
gene_annot0 <- gene_annot[, c("loc",
                              "loc_ID",
                              "chr",
                              "type",
                              "start",
                              "end",
                              "strand")]
##################################################
#head(motif_gene_matches, n = 5)
#Filter motifs for ranges (again? - just to be sure)
filtered_motif_df <- motif_gene_matches[motif_gene_matches$mstart <= 1500 |
                                          abs(
                                            motif_gene_matches$mstart-(
                                              motif_gene_matches$gene_end - motif_gene_matches$gene_start)+1) <= 1500, ]
##################################################
# Filter the data frame based on Filter_motif_orient variable
if (Filter_motif_orient == "forward") {
  filtered_motif_df <- filtered_motif_df %>%
    filter(grepl("F_", motif))
} else if (Filter_motif_orient == "reverse") {
  filtered_motif_df0 <- filtered_motif_df %>%
    filter(grepl("R_", motif))
} else if (Filter_motif_orient == "none") {
  filtered_motif_df <- filtered_motif_df %>%
    filter(grepl("F_|R_", motif))
  # No filter applied
} else {
  # Handle cases where Filter_motif_orient has unexpected values
  stop("Invalid value for Filter_motif_orient")
}
##################################################
#nothing changed in the testset (good?)
subset_filmotif_df <- filtered_motif_df[, c("loc",
                                            "motif",
                                            "mstart",
                                            "mend",
                                            "score",
                                            "strand")]
#gene_annot0$loc <- gsub("^.*:", "", gene_annot0$loc)
#subset_filmotif_df$loc <- gsub("^.*:", "", subset_filmotif_df$loc)
#message("Description of chromosomes (e.g., versions) is reduced to numerals here.")


merg_df <- merge(gene_annot0,
                 subset_filmotif_df,
                 by = "loc", all.x = TRUE)
merg_df <- na.omit(merg_df)

# Internal check for the number of rows
if (nrow(merg_df) < 2) {
  error_message <- "Please check identifiers chr#:start-end for fasta and gff input"
  stop(error_message)
}

################################################## RC RC RC RC RC
##### mstart is already the distance to the gen position 1 and end!!!!!
##### Calculate the genomic location here
merg_df$gen_mstart <- merg_df$mstart -1000 + merg_df$start
# find distances to start and end of the region of interest
#merg_df$flank_end <- merg_df$end + 1000
merg_df$dist_gen_start <- abs(
  merg_df$gen_mstart - merg_df$start +1000) 
merg_df$dist_gen_end <- abs(
  merg_df$end + 1000 - merg_df$gen_mstart )
# Create region column based on flank_start and flank_end
merg_df$region <- ifelse(
  merg_df$dist_gen_start  < 
    merg_df$dist_gen_end,
  "upstream",
  "downstream")
#############################################################################
#merg_df<- subset(merg_df, score >= "20")
############################################################################# THIS STEP IS NOT TRIVIAL
########################### TRANSCRIPTS THAT ARE SMALLER THAN 1K CAN HAVE MATCHES IN UP AND DOWN STREAM
############# REGIONS. THERE IS JUST ONE MATCH BUT IT IS EQUALLY REPRESENTED IN BOTH DATA SETS.
#### FILTER BY MOTIF PREFERENCES ########## MAY INCLUDE FLAG THAT MATCH APPEARS IN UP AND DOWN - LATER 
# Subsetting columns from merg_df to create merg_df1 without dist_gen_start ###########################
merg_df1 <- merg_df[, !(colnames(merg_df) %in% c("dist_gen_start"))]
merg_df2 <- merg_df[, !(colnames(merg_df) %in% c("dist_gen_end"))]
#merg_df1 <- merg_df1[merg_df$dist_gen_start <= 1500, ]
#merg_df2 <- merg_df2[merg_df$dist_gen_end <= 1500, ]
# Renaming columns dist_gen_start and dist_gen_end to dist_transc_border in merg_df1
colnames(merg_df1)[colnames(merg_df1) == "dist_gen_end"] <- "dist_transc_border"
colnames(merg_df1)[colnames(merg_df1) == "dist_gen_start"] <- "dist_transc_border"
# Renaming columns dist_gen_start and dist_gen_end to dist_transc_border in merg_df2
colnames(merg_df2)[colnames(merg_df2) == "dist_gen_start"] <- "dist_transc_border"
colnames(merg_df2)[colnames(merg_df2) == "dist_gen_end"] <- "dist_transc_border"
# Joining merg_df1 and merg_df2 back together to create merg_df00
merg_df01 <- rbind(merg_df1, merg_df2)
merg_df00 <- merg_df01[merg_df01$dist_transc_border <= 1500, ]
#############################################################################
# Invert region based on strand.x column
merg_df00$region <- ifelse(
  merg_df00$strand.x == "-",
  ifelse(
    merg_df00$region == "upstream",
    "downstream",
    "upstream"),
  merg_df00$region)
##### CREATE ONE OUTPUT FILE THAT CONTAINS UNSPECIFIC MATCHES ###############
#############################################################################
#############################################################################
##### START TO FILTER USING THE MOTIF SEQLET SUMMARY STATISTICS #############
merg_df00$epm <- convert_epm(merg_df00$motif)
##
## #### ####### #################### ######################## NOW COMPLICATED STUFF
up_merg_df00 <- merg_df00[merg_df00$region == "upstream", ]
#subset_up_merg <- up_merg_df00[, c(1, 14, 15, 16)]
subset_tss_m <- tss_motifs[, c(1, 2, 3, 7, 8)]
merg_up_mo <- merge(up_merg_df00,
                    subset_tss_m,
                    by = "epm", all.x = TRUE)
merg_up_mo0 <- na.omit(merg_up_mo)
# Internal check for the number of rows
if (nrow(merg_up_mo0) < 2) {
  error_message <- "Please check identifiers in epm for TSS/TTS ranges.csv and occurences.txt"
  stop(error_message)
}


merg_up_mima_mo0 <- merg_up_mo0 %>%
  filter(dist_transc_border >= min & dist_transc_border <= max)
merg_up_q10_q90_mo0 <- merg_up_mo0 %>% 
  filter(dist_transc_border >= q10 & dist_transc_border <= q90)
do_merg_df00 <- merg_df00[merg_df00$region == "downstream", ]
subset_tts_m <- tts_motifs[, c(1, 2, 3, 7, 8)]
subset_tts_m$min <- subset_tts_m$min -1520
subset_tts_m$max <- subset_tts_m$max -1520
subset_tts_m$q10 <- subset_tts_m$q10 -1520
subset_tts_m$q90 <- subset_tts_m$q90 -1520
merg_do_mo <- merge(do_merg_df00,
                    subset_tts_m,
                    by = "epm", all.x = TRUE)
merg_do_mo0 <- na.omit(merg_do_mo)
# Internal check for the number of rows
if (nrow(merg_do_mo0) < 2) {
  error_message <- "Please check identifiers in epm for TSS/TTS ranges.csv and occurences.txt"
  stop(error_message)
}

merg_do_mima_mo0 <- merg_do_mo0 %>%
  filter(dist_transc_border >= min & dist_transc_border <= max)
merg_do_q10_q90_mo0 <- merg_do_mo0 %>% 
  filter(dist_transc_border >= q10 & dist_transc_border <= q90)
a_mima_df01 <- rbind(merg_up_mima_mo0,
                     merg_do_mima_mo0)
a_q10q90_df01 <- rbind(merg_up_q10_q90_mo0,
                       merg_do_q10_q90_mo0)
a0_DF_df <- merg_df00[, c(2, 3, 5, 6, 7,8,15,16,9,14,13,1,11,12)]
a0_mima_df <- a_mima_df01[, c(2, 3, 5, 6, 7,8,15,16,9,14,13,1,11,12)]
a0_q1q9_df <- a_q10q90_df01[, c(2, 3, 5, 6, 7,8,15,16,9,14,13,1,11,12)]
##############################################################################
################################################## (similar to mo_finder line 139, expect for feature filter)
#Create loc to idx df
gene_annotX$loc <- paste(
  gene_annotX$chr,
  ":",
  gene_annotX$start - 1000,
  "-",
  gene_annotX$end + 1000,
  sep = "")
#Reduce attributes
gene_annotX$loc_ID <- str_extract(
  gene_annotX$attributes,
  "(?<=ID=).*?(?=;)")
gene_annotX$loc_ID <- ifelse(
  grepl(":",
        gene_annotX$loc_ID),
  sub(".*:", "",
      gene_annotX$loc_ID),
  gene_annotX$loc_ID)
#subset the df
gene_annotX <- gene_annotX[, c("loc",
                               "loc_ID",
                               "chr",
                               "type",
                               "start",
                               "end",
                               "strand")]

x <- gene_annotX[gene_annotX$type != "gene", ]
x <- x[, c("loc_ID",
           "loc",
           "type",
           "start",
           "end")]
x$loc_ID <- sub("\\..*", "", x$loc_ID)
y <- gene_annotX[gene_annotX$type == "gene", ]
############################################################################
#CREATE THE FEATURE "untranscribed"
gene_annotU <- gene_annot0
gene_annotU$type.a <- c("untranscr")
gene_annotU$start.a <- gene_annotU$start -1000
gene_annotU$end.a <- gene_annotU$start -1
gene_annotU$start.b <- gene_annotU$end +1
gene_annotU$end.b <- gene_annotU$end +1000
gene_annotU1 <- gene_annotU[, c("loc",
                                "loc_ID",
                                "chr",
                                "type.a",
                                "start.a",
                                "end.a",
                                "strand")]
colnames(gene_annotU1) <- c("loc",
                            "loc_ID",
                            "chr",
                            "type",
                            "start",
                            "end",
                            "strand")
gene_annotU2 <- gene_annotU[, c("loc",
                                "loc_ID",
                                "chr",
                                "type.a",
                                "start.b",
                                "end.b",
                                "strand")]
colnames(gene_annotU2) <- c("loc",
                            "loc_ID",
                            "chr",
                            "type",
                            "start",
                            "end",
                            "strand")

gene_annotU12 <- rbind(gene_annotU1, gene_annotU2)
z <- gene_annotU12[, c("loc_ID",
                       "loc",
                       "type",
                       "start",
                       "end")]
x0 <- rbind(x,z)

x0$loc_ID <- gsub("[-']+", "", x0$loc_ID)
x0$loc_ID <- gsub("\\..*", "", x0$loc_ID)
y$loc_ID <- gsub("[-']+", "", y$loc_ID)
y$loc_ID <- gsub("\\..*", "", y$loc_ID)
############################################################################
feat_annot0<- merge(y, x0, by = "loc")
feat_annot1 <- feat_annot0
feat_annot2 <- feat_annot1[, c("loc_ID.x",
                               "loc",
                               "type.y",
                               "start.y",
                               "end.y")]

colnames(feat_annot2) <- c("loc_ID",
                           "loc",
                           "type.y",
                           "start.y",
                           "end.y")

feat_annot2 <- na.omit(feat_annot2)
feat_annot2 <- feat_annot2[!duplicated(feat_annot2),]
a_mima_df01 <- na.omit(a_mima_df01)
a_mima_df01 <- a_mima_df01[!duplicated(a_mima_df01),]
##########################################################
a_mima_df01$loc_ID <- gsub("[-']+", "", a_mima_df01$loc_ID)
a_mima_df01$loc_ID <- gsub("\\..*", "", a_mima_df01$loc_ID)


feat_anno_mima <- merge (a_mima_df01,
                         feat_annot2,by = "loc")
feat_anno_mima0 <- na.omit(feat_anno_mima)
feat_anno_mima0 <- feat_anno_mima0[!duplicated(feat_anno_mima0),]

feat_anno_mima01 <- subset(feat_anno_mima0, gen_mstart > start.y & gen_mstart < end.y)
colnames(feat_anno_mima01)
feat_anno_mima01 <- feat_anno_mima01[, c("loc",
                                         "loc_ID.x",
                                         "type",
                                         "start",
                                         "end",
                                         "motif",
                                         "gen_mstart",
                                         "dist_transc_border",
                                         "region",
                                         "type.y",
                                         "min",
                                         "max",
                                         "q10",
                                         "q90")]
feat_anno_mima01 <- feat_anno_mima01[!duplicated(feat_anno_mima01),]
table(feat_anno_mima01$type.y)
############################################################## CONTINUE HERE
colnames(a0_DF_df)
a0_DF_df$gen_mend <- c(a0_DF_df$gen_mstart + 14)
a0_mima_df$gen_mend <- c(a0_mima_df$gen_mstart + 14)
a0_q1q9_df$gen_mend <- c(a0_q1q9_df$gen_mstart + 14)
############################################################## PLACEHOLDER
#a0_DF_df$score <- c(100)
#a0_mima_df$score <- c(100)
#a0_q1q9_df$score <- c(100)
a0_mima_df$chr <- sub(":.*", "", a0_mima_df$loc)
a0_q1q9_df$chr <- sub(":.*", "", a0_q1q9_df$loc)
a0_DF_bed <- a0_DF_df[, c("chr",
                          "gen_mstart",
                          "gen_mend",
                          "motif",
                          "score",
                          "strand.y")]
a0_mima_bed <- a0_mima_df[, c("chr",
                              "gen_mstart",
                              "gen_mend",
                              "motif",
                              "score",
                              "strand.y")]
a0_q1q9_bed <- a0_q1q9_df[, c("chr",
                             "gen_mstart",
                             "gen_mend",
                             "motif",
                             "score",
                             "strand.y")]

colnames(a0_DF_bed) <- c("chrom",
                         "chromStart",
                         "chromEnd",
                         "name",
                         "score",
                         "strand")
colnames(a0_mima_bed) <- c("chrom",
                           "chromStart",
                           "chromEnd",
                           "name",
                           "score",
                           "strand")
colnames(a0_q1q9_bed) <- c("chrom",
                          "chromStart",
                          "chromEnd",
                          "name",
                          "score",
                          "strand")

a0_DF_bed$thickStart <- c(a0_DF_bed$chromStart)
a0_DF_bed$thickStart <- c(a0_DF_bed$chromEnd)
a0_DF_bed$itemRgb <- c("255,0,127")

a0_mima_bed$thickStart <- c(a0_mima_bed$chromStart)
a0_mima_bed$thickStart <- c(a0_mima_bed$chromEnd)
a0_mima_bed$itemRgb <- c("255,0,127")

a0_q1q9_bed$thickStart <- c(a0_q1q9_bed$chromStart)
a0_q1q9_bed$thickStart <- c(a0_q1q9_bed$chromEnd)
a0_q1q9_bed$itemRgb <- c("255,0,127")

############################################################## GENERATE OUTPUT
file_path_out <- file.path(dirpath_2, paste0(SPEC, MODEL, PROJECT,
                                             "_",Filter_annot_type,"_", Filter_motif_orient
                                             , DATA))
#write.csv(a0_DF_df, file=paste0(file_path_out,
#                                "-def.csv"), row.names=FALSE)
write.csv(a0_mima_df, file=paste0(file_path_out,
                                  "-mima.csv"), row.names=FALSE)
write.csv(a0_q1q9_df, file=paste0(file_path_out,
                                  "-q1q9.csv"), row.names=FALSE)

#write.table(a0_DF_bed, file=paste0(file_path_out,"-def.bed"), row.names=FALSE, sep = "\t", quote = FALSE)
write.table(a0_mima_bed, file = paste0(file_path_out, "-mima.bed"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(a0_q1q9_bed, file=paste0(file_path_out, "-q1q9.bed"), row.names=FALSE, sep = "\t", quote = FALSE)

write.csv(feat_anno_mima01, file=paste0(file_path_out,"feat_mima_q1q9.csv"), row.names = FALSE)
##############################################################
# Step 1: motif names in file 1, 2 (_motifs, [epm]: epm_Soly_M_p0m00) 
# must match to file 4 out_all (motif_gene_matches, [motif]: epmSola_S035_p1m04F_640_16.2_AAAAAAWAHAAAA)
# Step 2: Position of genes must be normalized, file 3 (gene_annot, [start]-1000, [end]+1000) 
# pasted in new col "loc" with [chr]
# must match to file 4 out all (motif_gene_matches, [gene_start], [gene_end])
# Step 3: The file 3 and file 4 can be joined (reduced to essential number of columns)
# Step 4: Determine the orientation of genes according to [strand] and 
# re-calculate the gene/motif start and gene/motif end positions. 
# Step 5: Filter motifs based on file 1 and 2