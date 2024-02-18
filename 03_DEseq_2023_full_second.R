library(tximport)
library(GenomicFeatures)
library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(gplots)

rm(list = ls())
dirs <- list.files("second_batch/Trinity_full_salmon_quant/", "_quant")
files <- file.path("second_batch/Trinity_full_salmon_quant", dirs, "quant.sf")
names(files) <- dirs 
######### gene level
tx2gene <- read_delim(file.path("second_batch/Trinity_full_salmon_quant/", "Trinity.fasta.gene_trans_map"),
                      col_names = FALSE) %>% 
    dplyr::select(X2, X1)
names(tx2gene) <- c("transcript_IDs", "gene_IDs")
txi.salmon.g <- tximport(files, type = "salmon", tx2gene=tx2gene)

# meta information
samples <- read.table("second_batch/meta.txt", header = TRUE)
samples$color <- factor(samples$color)
samples$color <- relevel(samples$color, ref = "Brown")
all(str_sub(dirs, end = 6) == str_sub(samples$individual, end = 6))

dds_g <- DESeqDataSetFromTximport(txi.salmon.g, samples, ~color)
keep <- rowSums(counts(dds_g) >= 10 ) >= 10 
dds_g <- dds_g[keep,]
dds_g <- DESeq(dds_g)
res_g <- results(dds_g, alpha = 0.05)
# > resultsNames(dds_g)
# [1] "Intercept"            "color_Green_vs_Brown"
summary(res_g)
plotDispEsts(dds_g, main="Dispersion plot")
plotMA(res_g, alpha = 0.05) # this alpha is FDR; not the ajusted p-value in the table

resLFC <- lfcShrink(dds_g, coef="color_Green_vs_Brown", type="apeglm")
plotMA(resLFC, alpha = 0.05)

png("second_batch_full_trinity_volcano.png", width = 2100, height = 2100, units = "px", res = 300)
par(mfrow=c(1,1))
with(res_g, plot(log2FoldChange, -log10(pvalue), pch=20, 
               xlim=c(-12,12)))
with(subset(res_g, padj<.05),
     points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res_g, padj<.05 & abs(log2FoldChange) > 2), 
     points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res_g, padj<.05 & log2FoldChange > 2), 
     points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
dev.off()

LFC2up <- res_g[res_g$log2FoldChange > 2,]
LFC2upOrdered <- LFC2up[order(desc(LFC2up$log2FoldChange), LFC2up$padj), ]
write.table(as.data.frame(LFC2upOrdered),
          file="second_batch/trinity_full_LFC2upOrdered_results.txt",
          sep = "\t")

up_filters <- which(res_g@listData$padj<0.05 & res_g@listData$log2FoldChange > 0)
up_df <- data.frame('trans'=res_g@rownames[up_filters], 
                    'log2FoldChange' = res_g@listData$log2FoldChange[up_filters],
                    'pvalue'=res_g@listData$padj[up_filters])
up_df

LFC2down <- res_g[res_g$log2FoldChange < -2,]
LFC2downOrdered <- LFC2down[order(LFC2down$log2FoldChange, LFC2down$padj), ]
write.table(as.data.frame(LFC2downOrdered),
            file="second_batch/trinity_full_LFC2downOrdered_results.txt",
            sep = "\t")

down_filters <- which(res_g@listData$padj<0.05 & res_g@listData$log2FoldChange < 0)
down_df <- data.frame('trans'=res_g@rownames[down_filters], 
                    'log2FoldChange' = res_g@listData$log2FoldChange[down_filters],
                    'pvalue'=res_g@listData$padj[down_filters])
down_df
