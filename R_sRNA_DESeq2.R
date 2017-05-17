library(DEseq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(reshape2)
directory <- "~/smallrna_timecourse/nonmirs/HTSeq-DESeq2/5kb_counts/" # Set this to dir with JUST counts files.
sampleFiles <- list.files(directory)
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, time = rep(c("0", "0", "1", "6", "12", "24"), times = 2))
de.e.table <- sampleTable[-c(2,8),]
de.0.table <- sampleTable[-c(1,7),]
de.0 <- DESeqDataSetFromHTSeqCount(sampleTable = de.0.table, directory = directory, design = ~ time)
de.e <- DESeqDataSetFromHTSeqCount(sampleTable = de.e.table, directory = directory, design = ~ time)
de.e$time <- factor(de.e$time, levels = c("0", "1", "6", "12", "24"))
de.0$time <- factor(de.0$time, levels = c("0", "1", "6", "12", "24"))
de.0.LRT <- DESeq(de.0, test = "LRT", reduced = ~ 1)
de.e.LRT <- DESeq(de.e, test = "LRT", reduced = ~ 1)
de.e <- DESeq(de.e)
de.0 <- DESeq(de.0)
res.0 <- results(de.0)
res.0.LRT <- results(de.0.LRT)
res.e.LRT <- results(de.e.LRT)
res.e <- results(de.e)
res.0 <- res.0[order(res.0$pvalue),]
res.e <- res.e[order(res.e$pvalue),]
res.e.LRT <- res.e.LRT[order(res.e.LRT$pvalue),]
res.0.LRT <- res.0.LRT[order(res.0.LRT$pvalue),]
summary(res.0)
summary(res.0.LRT)
summary(res.e)
summary(res.e.LRT)
res.0 <- res.0[res.0$baseMean > 0,]
res.e <- res.e[res.e$baseMean > 0,]
res.e.LRT <- res.e.LRT[res.e.LRT$baseMean > 0,]
res.0.LRT <- res.0.LRT[res.0.LRT$baseMean > 0,]

df.0 <- data.frame(site = rownames(res.0), log2FoldChange = res.0$log2FoldChange,
                   pvalue = res.0$pvalue, padj = res.0$padj)
df.0 <- mutate(df.0, sig=ifelse(df.0$padj<0.1, "padj<0.1", "n.s."))
ggplot(df.0, aes(x = log2FoldChange, y = -log10(pvalue))) + geom_point(aes(col = sig)) + 
  scale_colour_manual(values = sciPalette) + theme_bw() +
  #geom_text_repel(data=filter(df.0, pvalue<0.05), aes(label=site)) +
  labs(title = "Volcano plot",
       subtitle = "0 hr (-OHT) vs 1 - 24 hr (+OHT)")

df.e <- data.frame(site = rownames(res.e), log2FoldChange = res.e$log2FoldChange,
                   pvalue = res.e$pvalue, padj = res.e$padj)
df.e <- mutate(df.e, sig=ifelse(df.e$padj<0.1, "padj<0.1", "n.s."))
ggplot(df.e, aes(x = log2FoldChange, y = -log10(pvalue))) + geom_point(aes(col = sig)) + 
  scale_colour_manual(values = sciPalette) + theme_bw() + 
  geom_text_repel(data=filter(df.e, pvalue<0.05), aes(label=site)) +
  labs(title = "Volcano plot",
       subtitle = "EV vs 1 - 24 hr (+OHT)")


