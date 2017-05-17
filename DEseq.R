#__miR DE__

# Only care about known miRNAs, probably high confidence (use this first).

# Number of miRNAs mapping to the hairpins, get a counts table, use DEseq2/edgeR.

# htseq-count can do this. Need hairpin GFFs, then it is very simple:
#$ htseq-count [options] <alignment_file> <gff_file>
  
  # see: http://www-huber.embl.de/HTSeq/doc/count.html#count
  
  # Got gff3 from ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3
  #This has the first column as chrN rather than my usual N. Used sed to change this:
  
#  $ cat hsa.gff3 | sed 's/^chr//' > hsa_matureANDpri.gff3

#$ for s in {01..12}; do htseq-count --format=bam --type=miRNA_primary_transcript --idattr=Name --quiet srna$s.srtd.bam ~/indexes/hsa_mirnas_matureANDpri.gff3 > primarycounts_srna$s.counts; done

#$ for s in {01..12}; do htseq-count --format=bam --type=miRNA --idattr=Name --quiet srna$s.srtd.bam ~/indexes/hsa_mirnas_matureANDpri.gff3 > maturecounts_srna$s.counts; done

# Will use DEseq2 most likely to apply statistical models to the count tables.
# In R now.

directory <- ""
sampleFiles <- list.files(directory)
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = rep(c("EV", rep("AsiSI", times = 5)), times = 2), time = rep(c("0", "0", "1", "6", "12", "24"), times = 2))
asisi.table <- subset(sampleTable, sampleTable$condition != "EV")
# This just selects for the AsiSI containing samples, as EV complicates the analysis.
de.asisi <- DESeqDataSetFromHTSeqCount(sampleTable = asisi.table, directory = directory, design = ~ time)
de.asisi$time <- factor(de.asisi$time, levels = c("0", "1", "6", "12", "24"))
de.asisi <- DESeq(de.asisi)
res <- results(de.asisi)
summary(res) # Get summary of results; tells us no. sig. changers up and down etc.
plotMA(res, ylim = c(-1.5,1.5)) # Get an MA plot for thesis fodder! Can add main etc.
resOrdered <- res[order(res$padj),]
plotCounts(de.asisi, gene = which.min(res$padj), intgroup = "time") # Check counts for most sig.
resSig <- subset(resOrdered, padj < 0.1) # Subset for only significant hits
write.csv(as.data.frame(resSig), file = "timecourseSIG.csv")

# Can use ggplot to create nicer graphs/boxplots. First use plotCounts to search for your miRNA of interest, then add the argument "returnData = T" and assign to a dataframe (e.g. sub.d):

ggplot(sub.d, aes(x = condition, y = count)) + geom_boxplot() + ggtitle("miR-34a") + labs(x = "Time post-OHT induction (hours)", y = "Normalized counts")

# Worked out how to create boxplots over a continous x axis for nice curves
# Need to ensure that the dataframe column for time is numeric

ggplot(new, aes(time, count)) + geom_boxplot(aes(group = time)) + stat_smooth(se = F, colour = "lightgrey", span = 1, size = 0.5, linetype = "dashed") + scale_x_continuous(breaks = c(0,1,6,12,24), "Time post-OHT induction (hours)") + scale_y_continuous("Normalized counts") + ggtitle("miRNA-34a expression") + theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA, colour = "black"))

