## Count matrix prepartation 
rm(list=ls())
setwd("/Volumes/Seagate_Backup_Plus_Drive/zika")
sampleTable <- read.table("sample_table.txt", header = T, sep = "\t")
filenames <- file.path(paste0(sampleTable$Run_s, ".bam"))
library("Rsamtools")
bamfiles <- BamFileList(filenames)
library("GenomicFeatures")
library("GenomicAlignments")
gtffile <- file.path("Mus_musculus.GRCm38.86.chr.gtf")
txdb <- makeTxDbFromGFF(gtffile, format="gtf")
genes <- exonsBy(txdb, by="gene")
genes <-renameSeqlevels(genes,c("1"="chr1", "2"="chr2", "3"="chr3", 
                                "4"="chr4","5"="chr5", "6"="chr6", 
                                "7"="chr7", "8"="chr8", "9"="chr9", 
                                "10"="chr10", "11"="chr11", 
                                "12"="chr12", "13"="chr13", 
                                "14"="chr14","15"="chr15", 
                                "16"="chr16", "17"="chr17", 
                                "18"="chr18","19"="chr19", 
                                "X"="chrX", "Y"="chrY", "MT"="chrM"))
se <- summarizeOverlaps(features=genes, reads=bamfiles, singleEnd=TRUE)
save(file="zika_objects_all.Rda", se)
load("zika_objects_all1.Rda")
head(assay(se))
se$infection_status_s
#check the millions of fragments that uniquely aligned to the genes. 
round(colSums(assay(se)) / 1e6, 1 )
## Pre-filtering the dataset and rlog transformation 
library(DESeq2)
dds <- DESeqDataSet(se, design = ~ infection_status_s)
nrow(dds)
#Pre-filtering the dataset
dds <- dds[rowSums(counts(dds)) > 1, ]
nrow(dds)
dds <- estimateSizeFactors(dds)
#transform to log2 scale
#sequencing depth correction is done automatically for the rlog method
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
#Scatterplot of transformed counts from two samples (Mock1 vs Zika1)
pdf("Scatterplot.pdf", width=6, height=3)
par( mfrow = c( 1, 2 ) )
plot(log2(counts(dds, normalized=TRUE)[,c(1,3)] + 1),
     pch=16, cex=0.3, xlab="not infected", ylab="ZIKV infected", main="log2 transform")
plot(assay(rld)[, c(1:3)],
     pch=16, cex=0.3, xlab="not infected", ylab="ZIKV infected", main="rlog transform")
dev.off()
#PCA plotting
library("pheatmap")
library("RColorBrewer")
sampleDists <- dist(t( assay(rld) ) )
sampleDistMatrix <- as.matrix(sampleDists )
rownames(sampleDistMatrix) <- paste( rld$infection_status_s, sep="-" )
pdf("PCA.pdf", width=6, height=4)
plotPCA(rld, intgroup = c("infection_status_s"))
dev.off()
#Differential expression analysis
##MA-plot
dds <- DESeq(dds)
pdf("MAplot.pdf", width=4.5, height=3)
plotMA(dds)
dev.off()
##Significant differentially expressed genes identified by DESeq2
res <- results(dds, alpha = 0.05)
mcols(res, use.names=TRUE)
summary(res)

#volcano plot
tab <- data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$padj))
pdf("Volcano.pdf", width=4.5, height=4)
par(mar = c(5, 4, 4, 5))
plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), 
     ylab = expression(-log[10]~adjpvalue), xlim=c(-5, 5), ylim=c(0, 60))
lfc = 1
pval = 0.05
signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
points(tab[signGenes, ], pch = 16, cex = 0.8, col = "red") 
abline(h = -log10(pval), col = "green3", lty = 2) 
abline(v = c(-lfc, lfc), col = "blue", lty = 2) 
mtext(paste("adjpval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
mtext(c(paste("-", lfc, "lfc"), paste("+", lfc, "lfc")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
dev.off()

#No. of significant genes
res.05 <- results(dds, alpha=.05)
table(res.05$padj <= .05)
#Fold change 
resLFC1 <- results(dds, lfcThreshold=1)
#Number of genes with Padj <= 0.05 and fold change >=2
table(resLFC1$padj <= 0.05)

##Significant differentially expressed genes identified by edgeR
library(edgeR)
edger <- DGEList(assay(se), group=colData(se)$infection_status_s)
edger <- edger[rowSums(cpm(edger)) > 1, ]
edger <- calcNormFactors(edger, method=c("TMM"))
edger <- estimateCommonDisp(edger)
edger <- estimateTagwiseDisp(edger)
res_edgeR <- exactTest(edger)
res_edgeR <- topTags(res_edgeR, n = nrow(edger$table), 
                     adjust.method="BH", p.value=1)$table
resSig_edgeR <- res_edgeR[res_edgeR$FDR < 0.05, ]
dim(resSig_edgeR)
##Method comparison: DESeq2 vs. edgeR
res_DESeq <- results(dds, alpha=.05)
res_DESeq <- res_DESeq[!is.na(res_DESeq$padj),]
resSig_DESeq <- res_DESeq[res_DESeq$padj < .05, ]
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R")
setlist <- list(edgeR=rownames(resSig_edgeR), DESeq=rownames(resSig_DESeq))
OLlist <- overLapper(setlist=setlist, sep="_", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)
pdf("Vennplot.pdf", width=4.5, height=4)
vennPlot(counts=counts, mysub="")
dev.off()

#Annotation
library("genefilter")
library("AnnotationDbi")
library("org.Mm.eg.db")
res$symbol <- mapIds(org.Mm.eg.db, keys=row.names(res), column="SYMBOL",
                     keytype="ENSEMBL", multiVals="first")
res$entrez <- mapIds(org.Mm.eg.db, keys=row.names(res), column="ENTREZID",
                     keytype="ENSEMBL",multiVals="first")
nres <- res[!is.na(res$symbol),]
resOrdered <- nres[order(nres$padj),]
resOrderedDF <- as.data.frame(resOrdered)
write.csv(resOrderedDF, file="results.csv")
sig <- resOrdered[!is.na(resOrdered$padj) &
                    resOrdered$padj<0.05 &
                    abs(resOrdered$log2FoldChange)>=1,]
write.table(rownames(sig), "sig.ensem.txt", row.names = FALSE, col.names=FALSE, quote=FALSE)
sig.symbol <- sig[, c(2:6)]
row.names(sig.symbol) <- sig$symbol
sig.symbol <- data.frame(sig.symbol)
write.table(rownames(sig.symbol), "sig.symbol.txt", row.names = FALSE, col.names=FALSE)
library(xtable)
print.xtable(xtable(sig.symbol, digits=c(0, 3,3, 3,-3, -3)), type = "latex", file = "siggenes.tex", 
             latex.environments = "center", include.rownames = TRUE)
#heatmapping
#select genes
selected <- rownames(sig)
selecteddds <- dds[rownames(dds) %in% selected,]
selectedrld <- rlog(selecteddds, blind=FALSE)
selectedmat <- assay(selectedrld)
selectedmat <- selectedmat - rowMeans(selectedmat)
selectedmat1 <- selectedmat
rownames(selectedmat1) <- sig$symbol[match(rownames(selectedmat1), rownames(sig))]
df <- as.data.frame(colData(selectedrld)[,c("infection_status_s")])
colnames(df) <- "group"
pdf("Heatmap.pdf", width=4.5, height=3, onefile=FALSE)
pheatmap(selectedmat1, annotation_col=df)
dev.off()
#Goterm enrichment analysis
universe <- rownames(resOrdered)
library(org.Mm.eg.db)
genemap <- select(org.Mm.eg.db, keys = selected,
                  columns=c("ENTREZID","SYMBOL","GENENAME"), keytype="ENSEMBL")
univmap <- select(org.Mm.eg.db, keys = universe,
                  columns=c("ENTREZID","SYMBOL","GENENAME"), keytype="ENSEMBL")
library(GOstats)
param<- new ("GOHyperGParams", geneIds = genemap, universeGeneIds=univmap, annotation="org.Mm.eg.db", ontology="BP",pvalueCutoff=0.01, conditional=FALSE, testDirection="over")
# run analysis
hyp<-hyperGTest(param)
# visualize
hypres <- summary(hyp)[, c(1, 2, 3, 7)]
print.xtable(xtable(hypres), type = "latex", file = "hypres.tex", 
             latex.environments = "center", include.rownames = FALSE)
