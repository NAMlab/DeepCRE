

# 0. Load packages ----

library(tidyverse)
library(edgeR) 
library(matrixStats) 
library(cowplot)
library(ggplot2)
library(plotly)
library(tools)




# 1. Import and organize data ----
# Required input files:
ssr_predictions_on_new_assemblies.csv # ssr model predictions for the 15 accessions
msr_predictions_on_new_assemblies.csv # msr model predictions for the 15 accessions
msr_predictions_on_pennellii.csv # msr model predictions for S. pennellii
itag_targets.csv # true labels for S. lycopersicum
spennellii_targets.csv # true labels for S. pennellii
mapman_lycopersicum_ITAG4.1.txt # Mercator4-generated MapMan annotation for S. lycopersicum
mapman_sopen.txt # Mercator4-generated MapMan annotation for S. pennellii

# Import the ssr data for the 15 accessions
dataset.ssr = read.csv("ssr_predictions_on_new_assemblies.csv", header = T, sep = ",", row.names = 1)
rownames(dataset.ssr) = gsub("\\..", "", rownames(dataset.ssr))

# Import the msr data for the 15 accessions
dataset.msr = read.csv("msr_predictions_on_new_assemblies.csv", header = T, sep = ",", row.names = 1)
rownames(dataset.msr) = gsub("\\..", "", rownames(dataset.msr))
dataset.msr.lyc = data.frame(lycopersicum = dataset.msr[,1], row.names = rownames(dataset.msr))

# Import the msr data for the pennellii
dataset.msr.pen = read.csv("msr_predictions_on_pennellii.csv", header = T, sep = ",", row.names = 1)
rownames(dataset.msr.pen) = gsub("\\..", "", rownames(dataset.msr.pen))


# Import the true label data for lycopersicum solgenomics
dataset.lyc = read.csv("itag_targets.csv", header = T, sep = ",", row.names = 1)
rownames(dataset.lyc) = gsub("\\..", "", rownames(dataset.lyc))

# Import the true label data for pennellii
dataset.pen = read.csv("spennellii_targets.csv", header = T, sep = ",", row.names = 1)
rownames(dataset.pen) = gsub("\\..", "", rownames(dataset.pen))

# 2. Compare between labels and predictions for the mrs models on lyc and pen ----

# Comparison of the predictions and true labels for lycopersicum
tmp = intersect(rownames(dataset.msr.lyc), rownames(dataset.lyc))
boxplot(dataset.msr.lyc[tmp,"lycopersicum"]~dataset.lyc[tmp,"true_target"])


# Comparison of the predictions and true labels for pennellii
tmp = intersect(rownames(dataset.msr.pen), rownames(dataset.pen))
boxplot(dataset.msr.pen[tmp,"pennellii"]~dataset.pen[tmp,"true_target"])

tmp = intersect(rownames(dataset.msr.lyc), rownames(dataset.lyc))
dataTibble.lyc = tibble(geneID = tmp, Prediction = dataset.msr.lyc[tmp,"lycopersicum"], TrueLabel = dataset.lyc[tmp,"true_target"])
dataTibble.lyc$TrueCode = c(0, 1, 0.5)[dataTibble.lyc$TrueLabel+1]
dataTibble.lyc$TrueLabel = factor(c("low", "high", "middle")[dataTibble.lyc$TrueLabel+1], levels= c("low", "middle", "high"))

lyc.TP = sum(dataTibble.lyc$TrueCode == 1 & dataTibble.lyc$Prediction > 0.5)
lyc.TN = sum(dataTibble.lyc$TrueCode == 0 & dataTibble.lyc$Prediction < 0.5)
lyc.FP = sum(dataTibble.lyc$TrueCode == 0 & dataTibble.lyc$Prediction > 0.5)
lyc.FN = sum(dataTibble.lyc$TrueCode == 1 & dataTibble.lyc$Prediction < 0.5)

lyc.Acc = sum(lyc.TP, lyc.TN)/ sum(lyc.TP,lyc.TN,lyc.FP,lyc.FN)
lyc.TPR = lyc.TP/sum(lyc.TP, lyc.FP)
lyc.TNR = lyc.TN/sum(lyc.TN, lyc.FN)

tmp = intersect(rownames(dataset.msr.pen), rownames(dataset.pen))
dataTibble.pen = tibble(geneID = tmp, Prediction = dataset.msr.pen[tmp,"pennellii"], TrueLabel = dataset.pen[tmp,"true_target"])
dataTibble.pen$TrueCode = c(0, 1, 0.5)[dataTibble.pen$TrueLabel+1]
dataTibble.pen$TrueLabel = factor(c("low", "high", "middle")[dataTibble.pen$TrueLabel+1], levels= c("low", "middle", "high"))

pen.TP = sum(dataTibble.pen$TrueCode == 1 & dataTibble.pen$Prediction > 0.5)
pen.TN = sum(dataTibble.pen$TrueCode == 0 & dataTibble.pen$Prediction < 0.5)
pen.FP = sum(dataTibble.pen$TrueCode == 0 & dataTibble.pen$Prediction > 0.5)
pen.FN = sum(dataTibble.pen$TrueCode == 1 & dataTibble.pen$Prediction < 0.5)

pen.Acc = sum(pen.TP, pen.TN)/ sum(pen.TP,pen.TN,pen.FP,pen.FN)
pen.TPR = pen.TP/sum(pen.TP, pen.FP)
pen.TNR = pen.TN/sum(pen.TN, pen.FN)

p1 = ggplot(dataTibble.lyc) +
  aes(x=TrueLabel, y=Prediction, fill=TrueLabel) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  ylim(0,1) +
  scale_fill_manual(values = c("#1573ba", "lightgoldenrod1", "#de425b")) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="Prediction ", x = "True Label",
       title="MSR performance",
       subtitle=expression(italic("S. lycopersicum"))) +
  theme_bw()


p2 = ggplot(dataTibble.pen) +
  aes(x=TrueLabel, y=Prediction, fill=TrueLabel) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  ylim(0,1) +
  scale_fill_manual(values = c("#1573ba", "lightgoldenrod1", "#de425b")) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="Prediction ", x = "True Label",
       title="MSR performance",
       subtitle=expression(italic("S. pennellii"))) +
  theme_bw()


pdf("Figure X (Prediction violin plots lycopersicum and pennellii.pdf",5,5)
plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
dev.off()

# 3. Loading MapMan bins ----

# lycopersicum
mapman.lyc = read.csv("mapman_lycopersicum_ITAG4.1.txt",quote = "", header = T, sep = "\t")[,1:3]
mapman.lyc = apply(mapman.lyc, 2, function(x) gsub("'", '', x, perl = TRUE))
mapman.lyc = as_tibble(mapman.lyc)
mapman.lyc$IDENTIFIER = gsub("\\..", "", mapman.lyc$IDENTIFIER)
mapman.lyc$IDENTIFIER = sapply(mapman.lyc$IDENTIFIER, toTitleCase)
mapman.lyc$BINCODE = paste("bin_",mapman.lyc$BINCODE,".", sep="")

tmp = match(mapman.lyc$IDENTIFIER,dataTibble.lyc$geneID)
mapman.lyc = bind_cols(mapman.lyc,dataTibble.lyc[tmp,])



mapman.lyc.bins = unique(mapman.lyc$BINCODE)
mapman.lyc.bins = sapply(mapman.lyc.bins,
                         function(x)
                           mapman.lyc[grep(x,mapman.lyc$BINCODE, fixed=T),,drop=F],
                         simplify = F)

tmp.new = sort(unlist(lapply(mapman.lyc.bins, nrow)),decreasing=T)
tmp.original[1:10]

mapman.bins = mapman.lyc[!duplicated(mapman.lyc$BINCODE),1:2]
mapman.bins = setNames(mapman.bins$NAME, mapman.bins$BINCODE)


rm(list=ls(pattern="tmp"))

# pennellii
mapman.pen = read.csv("mapman_sopen.txt",quote = "", header = T, sep = "\t")[,1:3]
mapman.pen = apply(mapman.pen, 2, function(x) gsub("'", '', x, perl = TRUE))
mapman.pen = as_tibble(mapman.pen)
mapman.pen$IDENTIFIER = gsub("\\..", "", mapman.pen$IDENTIFIER)
mapman.pen$IDENTIFIER = sapply(mapman.pen$IDENTIFIER, toTitleCase)
mapman.pen$BINCODE = paste("bin_",mapman.pen$BINCODE,".", sep="")

tmp = match(mapman.pen$IDENTIFIER,dataTibble.pen$geneID)
mapman.pen = bind_cols(mapman.pen,dataTibble.pen[tmp,])

mapman.pen.bins = unique(mapman.pen$BINCODE)
mapman.pen.bins = sapply(mapman.pen.bins,
                         function(x)
                           mapman.pen[grep(x,mapman.pen$BINCODE, fixed = T),,drop=F],
                         simplify = F)

rm(list=ls(pattern="tmp"))

# 4. Filtering the bins according to mapping coverage  ----
bins.nona = sapply(names(mapman.pen.bins), 
                   function(x){
                     sum(!is.na(mapman.lyc.bins[[x]]$Prediction)) > 2 &
                       sum(!is.na(mapman.pen.bins[[x]]$Prediction)) > 2
                   })

bins.nonasum = sapply(names(mapman.pen.bins), 
                      function(x){
                        mean(sum(!is.na(mapman.lyc.bins[[x]]$Prediction)),
                             sum(!is.na(mapman.pen.bins[[x]]$Prediction)))
                      })



mapman.pen.bins.filtered = mapman.pen.bins[bins.nona]
mapman.lyc.bins.filtered = mapman.lyc.bins[bins.nona]

# 5. Perform wilcoxon test to compare expression estimated between lyc and pen ----
bins.stats = sapply(names(mapman.lyc.bins.filtered), 
                    function(x){
                      as.vector(wilcox.test(
                        mapman.lyc.bins.filtered[[x]]$Prediction, 
                        mapman.pen.bins.filtered[[x]]$Prediction, 
                        paired = F, 
                        alternative = "t",
                        exact = F))
                    })
bins.stats = apply(bins.stats,2,unlist)

bins.stats = t(as.data.frame(bins.stats))
bins.stats = as_tibble(bins.stats)
bins.stats$FDR = p.adjust(bins.stats$p.value , method = "BH")
bins.stats$BINCODE = names(mapman.lyc.bins.filtered)

# 6. Filter bins according to the p-value threshold ----
bins.sig = bins.stats$BINCODE[bins.stats$FDR < 0.05]


mapman.significant = sapply(bins.sig,
                            function(x){
                              c(
                                BINCODE = x,
                                NAME = mapman.bins[x],
                                Median_LycPred = median(mapman.lyc.bins.filtered[[x]]$Prediction,na.rm=T),
                                Median_LycTrue = mean(mapman.lyc.bins.filtered[[x]]$TrueCode,na.rm=T),
                                Median_PenPred = median(mapman.pen.bins.filtered[[x]]$Prediction,na.rm=T),
                                Median_PenTrue = mean(mapman.pen.bins.filtered[[x]]$TrueCode,na.rm=T),
                                Diff_Pred = median(mapman.pen.bins.filtered[[x]]$Prediction,na.rm=T)-median(mapman.lyc.bins.filtered[[x]]$Prediction,na.rm=T),
                                Diff_True = mean(mapman.pen.bins.filtered[[x]]$TrueCode,na.rm=T)-mean(mapman.lyc.bins.filtered[[x]]$TrueCode,na.rm=T),
                                p.value = bins.stats$p.value[bins.stats$BINCODE == x],
                                FDR = bins.stats$FDR[bins.stats$BINCODE == x]
                              )
                            })
mapman.significant = t(as.data.frame(mapman.significant))
mapman.significant = as_tibble(mapman.significant)
mapman.significant$BINSIZE = as.numeric(bins.nonasum[mapman.significant$BINCODE])
mapman.significant$Label = sub("^.+\\.", "", mapman.significant$NAME.bin_1.1.)
mapman.significant$FDR = as.numeric(mapman.significant$FDR)
mapman.significant$Diff_True = as.numeric(mapman.significant$Diff_True)
mapman.significant$Diff_Pred = as.numeric(mapman.significant$Diff_Pred)
mapman.significant$p.value = as.numeric(mapman.significant$p.value)


# 7. Generate tables and figures ----

# Save results as a table
write.csv(mapman.significant, file= "Supplemental Table 7 (significant mapman bins).csv")


# Generate figures
p3 = ggplot(mapman.significant) + 
  aes(x = Diff_True, y = Diff_Pred, size = BINSIZE, label = Label, alpha = -log10(FDR) ) +
  geom_smooth(formula = y ~ x, method=lm, color = "darkgrey") +
  geom_point(shape=16, color = "dodgerblue4") +
  scale_size(range = c(2, 10)) +
  scale_alpha(range = c(0.2, 1)) +
  #geom_hex(show.legend = T, bins=20) +
  # coord_cartesian(xlim = seq(-0.1, 0.5, by = 0.1), ylim = seq(-0.1, 0.5, by = 0.1)) +
  scale_x_continuous(breaks = seq(-0.1, 0.5, by = 0.1)) +
  scale_y_continuous(breaks = seq(-0.1, 0.5, by = 0.1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #geom_text(data = subset(mapman.significant, -log10(FDR) > 1), aes(label = Label), hjust = 0, vjust = 0.1, color = "red") +
  #scale_shape_manual(values = 16) +
  labs(y="Predicted expression difference", x = "Observed expression difference",
       title="Differential expression of gene ontology terms",
       subtitle=expression(italic("S. pennellii - S. lycopersicum")),
       caption="Based on automatic Mercator4 MapMan annotation") +
  theme_classic()
# theme_dark() 
#theme_bw()

p12 = plot_grid(p1,p2, labels = c('a', 'b'), label_size = 12)

pdf("Figure X (Prediction violin plots lycopersicum and pennellii) v3.pdf",5,8)
plot_grid(p12,p3 + theme(legend.position = "none"), 
          nrow = 2,
          labels = c('', 'c'),
          rel_heights = c(1,2),
          label_size = 12)
dev.off()

pdf("Figure X legend.pdf",5,5)
plot_grid(p3 , 
          label_size = 12)
dev.off()

# make an interactive version of the ggplot scatterplot with plotly
p31 = ggplotly(p3, tooltip = c("label", "alpha"))
p31 <- layout(p31, hovermode = "closest")
p31


