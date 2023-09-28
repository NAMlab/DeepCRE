#install.packages("BiocManager")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("DNABin")
#BiocManager::install("ggtree")
#install.packages("phangorn")
#library(TFBSTools)
#library(JASPAR2020)
library(universalmotif)
library(motifmatchr)
library(ape)
library(motifStack)
library(ade4)
#library(phangorn)
library(ggtree)
##################################################
###################### Setup for "moca_blue" enviroment
NAME0="rdf5_epm"
SPEC ="Zema"
MODEL="C0" # C0 stand for DeepCistrome version 1 (available at 02-may-2023) Standard conditions
TYPE ="_pfm-motifs.jaspar"
#######################################################
#FILE1 = paste0(NAME0,SPEC,MODEL,TYPE)

FILE1 = "all_motifs_20230505.jaspar"
#######################################################
dirpath_in = "../Mo_Nom/out/"
dirpath_out = "./out/"
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
File1 <- paste0(dirpath_in,FILE1)
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################################################
cwm1 <- read_jaspar(File1)
##################################################
# Loop through each motif object in the list
for (i in seq_along(cwm1)) {
  
  # Extract the motif name and the number after the last "_" underscore
  motif_name <- attr(cwm1[[i]], "name")
  nsites <- as.numeric(sub(".+_(\\d+)$", "\\1", motif_name))
  
  # Assign the nsites value to the "nsites" field of the motif object
  cwm1[[i]]["nsites"] <- nsites
}
##################################################
# Loop through each motif object in the list
for (i in seq_along(cwm1)) {
  
  # Extract the Total IC and the Consensus values from the motif object
  total_ic <- attr(cwm1[[i]], "icscore")
  total_ic_rounded <- round(total_ic, 1)
  
  consensus <- attr(cwm1[[i]], "consensus")
  
  # Combine the Total IC and the Consensus values separated by "_" to the motif name
  motif_name <- attr(cwm1[[i]], "name")
  new_motif_name <- paste0(motif_name, "_", total_ic_rounded, "_", consensus)
  
  # Assign the new motif name to the motif object
  attr(cwm1[[i]], "name") <- new_motif_name
}
##################################################
pwm_uni0<-convert_motifs(
  cwm1, class = "TFBSTools-PWMatrix")

pcm<-convert_motifs(
  cwm1, class = "motifStack-pcm")
##################################################
sum<-as.data.frame(summarise_motifs(pcm))
write.csv(sum, file = paste0(dirpath_out,FILE1,"summary.txt"))
##################################################
c_pcm<-clusterMotifs(pcm) ### !!! TIME TO GET A COFFEEE !!! ###
hc<- c_pcm
motifs<-pcm[hc$order]
##################################################
write.tree(as.phylo(c_pcm), file = paste0(dirpath_out,FILE1,"-SW.nwk"))

########################################################
#motifs<-pcm[hc$order]
#motifs <- lapply(motifs, pcm2pfm)
#d1o alignment
#compare_motif()
##########################################################
cwm1[[1]]
##########################################################
c<-compare_motifs(cwm1, method = "EUCL")
c0<-as.data.frame(c)
write.csv(c0, file = paste0(dirpath_out,FILE1,"_matrix-EUCL.txt"))
#assign scores from data.frame to branches for selection
tree <- motif_tree(cwm1, layout = "rectangular", db.scores = "scores", method = "EUCL")
write.tree(as.phylo(tree), file = paste0(dirpath_out,FILE1,"-EUCL.nwk"))
##########################################################

##########################################################
c<-compare_motifs(cwm1, method = "WEUCL")
c0<-as.data.frame(c)
#assign scores from data.frame to branches for selection
tree <- motif_tree(cwm1, layout = "rectangular", db.scores = "scores", method = "WEUCL")
write.tree(as.phylo(tree), file = paste0(dirpath_out,FILE1,"-WEUC2.nwk"))
############################################################################
c<-compare_motifs(cwm1, method = "PCC")
c0<-as.data.frame(c)
#assign scores from data.frame to branches for selection
#filter scores based on failed IC (low motif scores)
tree <- motif_tree(cwm1, layout = "rectangular", db.scores = "scores", method = "PCC")
write.tree(as.phylo(tree), file = paste0(dirpath_out,FILE1,"-PCC.nwk"))
############################################################################
c<-compare_motifs(cwm1, method = "ALLR_LL")
c0<-as.data.frame(c)
#assign scores from data.frame to branches for selection
tree <- motif_tree(cwm1, layout = "rectangular", db.scores = "scores", method = "ALLR_LL")
write.tree(as.phylo(tree), file = paste0(dirpath_out,FILE1,"-ALLR_LL.nwk"))
############################################################################

#motif_pvalue(cwm1)



#ggtree::ggtree(b, 
#               layout="circular") +  geom_tiplab2(size =2.5) + xlim(-15, 5)
#cwm2 <- read_jaspar(File2)
#plot(hclust(dist(a0)))
#a<-compare_motifs(cwm1, method = "WPCC")
#a0<-as.data.frame(a)
#plot(hclust(dist(a0)))
#b<-hclust(dist(a0))

#ggtree::ggtree(b, 
#               layout="circular") +  geom_tiplab2(size =2.5) + xlim(-15, 5)
#b0<-as.phylo(b)
#write.tree(b0, file = paste(JOB_NAME,".nex"))

######################################################################



#motif_tree(cwm1, 
#           label = "name",
#           size=.3,
#           layout = "circular",
#           method = "EUCL",
#           legend = F,
#           db.scores = cwm1)


