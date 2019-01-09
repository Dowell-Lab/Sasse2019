###DESeq2 sample script
###margaret.gruca@colorado.edu, adapted in part from Stephen Turner, @genetics_blog

###Required libraries -----------------------------------------------------------------------------------
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
install.packages("gplots")
install.packages("ggplot2")
install.packages("RColorBrewer")
install.packages("calibrate")

###Different ways of pre-processesing count data : use only the lines you need ---------------------------

data_10 <- read.delim("/scratch/Users/magr0763/Gerber/GRO-seq/multibamcov/remap/Gerber_counts.hg38.multicov.bed", sep="\t", header=FALSE)
countdata <- write.table(data, file = "/scratch/your_directory/mysample.txt", append = FALSE, sep = "\t")
countdata <- write.csv(data, file = "/scratch/your_directory/mysample.txt")

###-------------------------------------------------------------------------------------------------------

## Previously ran at command line something like this:
## featureCounts -a genes.gtf -o counts.txt -T 12 -t exon -g gene_id GSM*.sam
#countdata <- read.table("counts.txt", header=TRUE, row.names=1)

# Remove first five columns (chr, start, end, strand, length)
#countdata <- countdata[ ,6:ncol(countdata)]

# Remove .bam or .sam from filenames
#colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))

#Import countdata from mutlibamcov/HTSeq txt/csv
data_30 <- read.csv("/scratch/Users/magr0763/Gerber/GRO-seq/Data_Analysis/DESeq/hg38_remap/count_data_30min.csv", header=TRUE, row.names=1)

# Convert to matrix -- countdata should only contain counts as body, samples as column names, and gene id as row names
countdata <- as.matrix(data_10)
countdata <- na.omit(countdata)
head(countdata)
#countdata <- na.omit(countdata)

# Assign condition (rep("condition", number of replicates)) -- must be in the same order as listed in countdata file!!
(condition <- factor(c(rep("Dex10", 2), rep("Untreated10", 2), rep("DexTNFa10", 2), rep("TNFa10", 2))))

# Analysis with DESeq2 ------------------------------------------------------------------------------------

library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)

#QC -- Check sizeFactors, should be close to 1 (ideally) -- largely affected by depth (read counts)
sizeFactors(dds)

#Can run option "iterate" which adjusts calculations assuming an ideal size facotor of 1
#dds <- estimateSizeFactors(dds, type= "iterate")

#Run the following to remove rows with number of reads < 100
#dds <- dds[ rowSums(counts(dds)) > 50, ] 

# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld), col="skyblue", border="slateblue", main="Gerber Histogram")

#Examine the effect of log transformation (QC check)

library(calibrate)

plot( log2( 1 + counts(dds, normalized=TRUE)[ , 1:2] ),
      col=rgb(0,0,0,.2), pch=16, cex=0.3 )
plot( assay(rld)[ , 1:2],
      col=rgb(0,0,0,.2), pch=16, cex=0.3 )

# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- dist( t( assay(rld) ) )
library(gplots)
png("qc-heatmap-samples2.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# Principal components analysis
## Could do with built-in DESeq2 function:
## DESeq2::plotPCA(rld, intgroup="condition")
## I like mine better:
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20) }

  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))

png("qc-pca2.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
dev.off()

condx <- "DexTNFa30"
condy <- "TNFa30"

# Get differential expression results, specify which conditions to contrast (max=2)
res <- results(dds, contrast = c("condition", condx, condy))
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file="diffexpr-results_DESeq2.csv")
write.table(resdata, file="diffexpr-results_DESeq2.txt", append = FALSE, sep = "\t" )

#Separate upregulated and downregulated genes w/padj cutoff

res <- na.omit(res)
resSig <- res[ res$padj < .05, ]
up <- subset(resSig, log2FoldChange > 0)
up <- merge(as.data.frame(up), as.data.frame(counts(dds,normalized=TRUE)), by="row.names", sort=FALSE)
write.table(up[order(up$log2FoldChange, decreasing=TRUE),],
            "UpRegulated.txt")
write.csv(up[order(up$log2FoldChange, decreasing=TRUE),],
            "UpRegulated.csv")
down <- subset(resSig, log2FoldChange < 0)
down <- merge(as.data.frame(down), as.data.frame(counts(dds,normalized=TRUE)), by="row.names", sort=FALSE)
write.table(down[order(down$log2FoldChange, decreasing=TRUE),],
            "DownRegulated.txt")
write.csv(down[order(down$log2FoldChange, decreasing=TRUE),],
            "DownRegulated.csv")


## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

## Check independent filtering
metadata(res)$filterThreshold

## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)

png("DE-MAplot.DESeq2.png", 1500, 1000, pointsize=20)
DESeq2::plotMA(dds, ylim=c(-6,6), cex=1, main="MA Plot")
dev.off()

## Plot with more detail:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="TNFa vs. Dex+TNFa : 30min", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-5.5, 5.5), ylim=c(0,80))
dev.off()

## Volcano plot with "significant" genes not labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="TNFa vs. Dex+TNFa : 30min", legendpos="bottomright", labelsig=FALSE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot2.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-5.5, 5.5), ylim=c(0,80))
dev.off()