library(magrittr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(Rgraphviz)
library(org.At.tair.db)

read.counts <- read.table("stringtie_development_transcript_count.csv", header = TRUE, sep=",")
row.names(read.counts) <- read.counts$transcript_id
head(read.counts)
read.counts <- read.counts[,-c(1)] #exclude all information that do not include counts
names(read.counts) <- c("GSM2879535",
                        "GSM2879536",
                        "GSM2879537",
                        "GSM2879538",
                        "GSM2879539",
                        "GSM2879540",
                        "GSM2879541",
                        "GSM2879542",
                        "GSM2879543",
                        "GSM2879553",
                        "GSM2879554",
                        "GSM2879555",
                        "GSM2879565",
                        "GSM2879566",
                        "GSM2879567",
                        "GSM2879577",
                        "GSM2879578",
                        "GSM2879579")
head(read.counts)

growscale_contol.groups <- c("100","100","100",
                            "200","200","200",
                            "300","300","300",
                            "400","400","400",
                            "500","500","500",
                            "600","600","600")

sample_info <- data.frame(condition=growscale_contol.groups,
                          row.names=names(read.counts))
DESeq.ds <- DESeqDataSetFromMatrix(countData = read.counts,
                                   colData = sample_info,
                                   design = ~ condition)
#colData(DESeq.ds)
#assay(DESeq.ds) %>% head
DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds))>0,]#eleminate counts that are less than zero
length(rowRanges(DESeq.ds)) # number of genes counted
DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds))>100,] #eleminate counts that are less than 200
names(rowRanges(DESeq.ds))
head(assay(DESeq.ds))
library(DESeq2)
reorder
#normalization
DESeq.ds <- estimateSizeFactors(DESeq.ds)
sizeFactors(DESeq.ds)
counts.sf_norm <- counts(DESeq.ds, normalized = TRUE)
#transformation
log.norm.counts <- log2(counts.sf_norm+1)
par(mar=c(9,5,3,3))
boxplot(counts.sf_norm, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(log.norm.counts , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (log.norm.counts [ ,8:9] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# r.log transformation ##########
library(vsn)
DESeq.rlog <- vst(DESeq.ds, blind = TRUE)
rlog.norm.counts <- assay(DESeq.rlog)
msd_plot <- meanSdPlot (log.norm.counts,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(rlog.norm.counts,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

# hierarchical clustering
distance.m_rlog <- as.dist(1- cor(rlog.norm.counts, method ="pearson"))
tiff("HC_makale.tiff", units="in", width=5, height=5, res=1200)
plot(hclust(distance.m_rlog),
      labels = colnames(rlog.norm.counts),
      xlab = "",
      main = "Hierarchical Clustering")
dev.off()
########### PCA ########################################
pc <- prcomp (t(rlog.norm.counts))
tiff("PCA_makale.tiff", units="in", width=5, height=5, res=1200)
plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(DESeq.ds)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(DESeq.rlog)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis")

print (P)
dev.off()
############## DEG analysis ############################
DESeq.ds <- DESeq(DESeq.ds)
DESeq.ds <- estimateSizeFactors(DESeq.ds)
DESeq.ds <- estimateDispersions(DESeq.ds)
DESeq.ds <- nbinomWaldTest(DESeq.ds)

DESeq.ds.results <- results(DESeq.ds, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 1)

write.table(as.data.frame(DESeq.ds.results), file = "DESeq.ds.results.txt",sep = "\t", quote = FALSE, row.names = TRUE)

DEG.names <- rownames(subset(DESeq.ds.results, 
                             padj < 0.05 & abs(log2FoldChange) > 1))
length(DEG.names)

DESeq.ds.results_200 <- results(DESeq.ds, 
                            independentFiltering = TRUE, 
                            alpha = 0.05,
                            contrast=c("condition","200","100"), lfcThreshold = 1)
#contrast=c("condition","300","100")
DESeq.ds.results_300 <- results(DESeq.ds, 
                            independentFiltering = TRUE, 
                            alpha = 0.05,
                            contrast=c("condition","300","100"), lfcThreshold = 1)
DESeq.ds.results_400 <- results(DESeq.ds, 
                                independentFiltering = TRUE, 
                                alpha = 0.05,
                                contrast=c("condition","400","100"),lfcThreshold = 1)
DESeq.ds.results_500 <- results(DESeq.ds, 
                                independentFiltering = TRUE, 
                                alpha = 0.05,
                                contrast=c("condition","500","100"), lfcThreshold = 1)
DESeq.ds.results_600 <- results(DESeq.ds, 
                                independentFiltering = TRUE, 
                                alpha = 0.05,
                                contrast=c("condition","600","100"), lfcThreshold = 1)


#resultsNames(DESeq.ds)
summary(DESeq.ds.results_600)
head(DESeq.ds.results$log2FoldChange)
global <- rownames(subset(DESeq.ds.results, padj < 0.05 & abs(log2FoldChange) > 1 ))
length(global)
one2two <- rownames(subset(DESeq.ds.results_200, padj < 0.05, abs(log2FoldChange) > 1))
length(one2two)
head(one2two)
one2three <- rownames(subset(DESeq.ds.results_300, padj < 0.05 & abs(log2FoldChange) > 1))
length(one2three)
one2four <- rownames(subset(DESeq.ds.results_400, padj < 0.05 & abs(log2FoldChange) > 1))
length(one2four)
one2five <- rownames(subset(DESeq.ds.results_500, padj < 0.05 & abs(log2FoldChange) > 1))
length(one2five)
one2six <- rownames(subset(DESeq.ds.results_600, padj < 0.05 & abs(log2FoldChange) > 1))
length(one2six)
library(VennDiagram)
venn.diagram(
  x = list(one2two, one2three, one2four, one2five, one2six),
  category.names = c("Set 1" , "Set 2 " , "Set 3", "Set 4", "Set 5"),
  filename = '#14_venn_diagramm.png',
  output=TRUE
)
tiff("Venn_diagram_makale.tiff", units="in", width=5, height=5, res=1200)
venn.plot <- venn.diagram(
  x = list(
    A = one2two,
    B = one2three,
    C = one2four,
    D = one2five,
    E = one2six
  ),
  filename = "Venn_5set-fold_REAL.tiff",
  col = "black",
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  alpha = 0.50,
  cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
          1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
  cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.cex = 1.5,
  cat.fontface = "bold",
  margin = 0.05
)
dev.off()
library(dplyr)
find.uniques <- list(a=one2two,
                       b=one2three,
                       c=one2four,
                       d=one2five,
                       e=one2six)
as <- find.uniques$a[!(find.uniques$a %in% c(find.uniques$b,
                                   find.uniques$c,
                                   find.uniques$d,
                                   find.uniques$e))]
bs <- find.uniques$b[!(find.uniques$b %in% c(find.uniques$a,
                                             find.uniques$c,
                                             find.uniques$d,
                                             find.uniques$e))]
cs <- find.uniques$c[!(find.uniques$c %in% c(find.uniques$a,
                                             find.uniques$b,
                                             find.uniques$d,
                                             find.uniques$e))]
ds <- find.uniques$d[!(find.uniques$d %in% c(find.uniques$a,
                                             find.uniques$b,
                                             find.uniques$c,
                                             find.uniques$e))]
es <- find.uniques$e[!(find.uniques$e %in% c(find.uniques$a,
                                             find.uniques$b,
                                             find.uniques$c,
                                             find.uniques$d))]


specific.genes <- as+bs+cs+ds+es
specific.genes
length(specific.genes)

DEG.names_300 <- (rownames(subset(DESeq.ds.results_300, padj < 0.05)))
length(DEG.names_300)

