library(DiffBind) #v3.8.1

# Save results and plots to PDF 
pdf("SN4741_37v39.pdf")

# Make the sample sheet containing all 37, 39, MB, and FB libraries
SampleID <- c("37-1", "37-2", "37-3", "37-4", "39-1", "39-2", "39-3", "39-4")
Condition <- c(rep("t37", 4), rep("t39", 4)) # add character in front of 37 and 39 because DiffBind doesn't love mask names that are numbers not characters in the future steps
Replicate <- c(rep(seq(1:4), 2))
bamReads <- c("37-1_SN4741_S1.bam", "37-2_SN4741_S3.bam", "37-3_SN4741_S2.bam", "37-4_SN4741_S8.bam", "39-1_SN4741_S4.bam", "39-2_SN4741_S7.bam", "39-3_SN4741_S5.bam", "39-4_SN4741_S6.bam")
Peaks <- c("37-1_SN4741_S1_peaks_blacklisted.bed", "37-2_SN4741_S3_peaks_blacklisted.bed", "37-3_SN4741_S2_peaks_blacklisted.bed", "37-4_SN4741_S8_peaks_blacklisted.bed", "39-1_SN4741_S4_peaks_blacklisted.bed", "39-2_SN4741_S7_peaks_blacklisted.bed", "39-3_SN4741_S5_peaks_blacklisted.bed", "39-4_SN4741_S6_peaks_blacklisted.bed")
PeakCaller <- c(rep("bed", 8))
sampleData_SN4741 <- data.frame(cbind(SampleID, Condition, Replicate, bamReads, Peaks, PeakCaller))

## PCA and correlation heatmaps of different stages of data analysis

# Load the data into Diffbind
# Consensus peakset: peak must be present in ≥2 libraries
SN4741_peaks <- dba(sampleSheet=sampleData_SN4741, minOverlap = 2)

# No read counts - just peak presence/absence
dba.plotHeatmap(SN4741_peaks) 
dba.plotPCA(SN4741_peaks, DBA_CONDITION, label = DBA_ID)

# Count reads overlapping peaks and do not normalize
SN4741_counts <- dba.count(SN4741_peaks, summits = 100, bRemoveDuplicates = TRUE, score=DBA_SCORE_READS)
save.image() # count takes a long time, save.

# Plot heatmap and PCA considering raw read counts
dba.plotHeatmap(SN4741_counts, score = DBA_SCORE_READS)
dba.plotPCA(SN4741_counts, score = DBA_SCORE_READS, DBA_CONDITION, label = DBA_ID)

# Normalize the read counts
# library size, RLE normalization
SN4741_normalized <- dba.normalize(SN4741_counts, library = "full", normalize = "RLE",  method=DBA_DESEQ2) 

# Tried background normalization - was a lot more skew in the MA plot that follows with background than full, particularly at those peaks with low read counts. Not sure why this difference exists since the two settings are meant to be proxies for each other but move forward with this method. Also go back to the MB and SN4741 comparison and change normalization to match. Since there were no explicit comparisons, didn't look at the MA plot up there, but expect the normalization to have behaved similiarly. Keep it consistent.

# Plot heatmap and PCA considering normalized read counts
dba.plotHeatmap(SN4741_normalized, score = DBA_SCORE_NORMALIZED)
dba.plotPCA(SN4741_normalized, score = DBA_SCORE_NORMALIZED, DBA_CONDITION, label = DBA_ID)

## Venn diagrams 

# How well do the replicates overlap? 
dba.plotVenn(SN4741_peaks, SN4741_peaks$masks$t37)
dba.plotVenn(SN4741_peaks, SN4741_peaks$masks$t39)

# How much do the two conditions overlap? 
venn <- dba.peakset(SN4741_peaks, consensus=DBA_CONDITION, minOverlap=2) 
dba.plotVenn(venn, venn$masks$Consensus)

## Differential accessbility + volcano plots 

# Set contrasts
SN4741_normalized <- dba.contrast(SN4741_normalized, contrast=c("Condition","t39","t37"))

# Analyse
SN4741_analyzed <- dba.analyze(SN4741_normalized, method = DBA_DESEQ2)
dba.analyze(SN4741_analyzed, bRetrieveAnalysis=TRUE)

# Check MA plot for weird distriubtions for normalization issues
# MA plots 
## all peaks, not normalized
dba.plotMA(SN4741_counts, contrast=list(Differentiated=SN4741_counts$masks$t39), bNormalized=FALSE,
sub="Non-Normalized")
## Normalized
dba.plotMA(SN4741_analyzed,method=DBA_DESEQ2)
dba.plotMA(SN4741_analyzed,method=DBA_DESEQ2, bSignificant = FALSE) # no pink overlay so can see distrib better

# Volcano plot
dba.plotVolcano(SN4741_analyzed, method=DBA_DESEQ2)

# Remaking the volcano plot so can customize the plot

## Write out results of analysis in way that's readable for native R plotting (or ggplot if you're feeling fancy)
SN4741_temperature.report <- dba.report(SN4741_analyzed, th = 1, fold = 0, method = DBA_DESEQ2, bCounts = T) # write out all peaks
SN4741_temperature.df <- as.data.frame(SN4741_temperature.report)

# Write out as a file
write.table(SN4741_temperature.df, "37v39_DAR_table.txt", quote = F)

## Replicating the volcano plot
with(SN4741_temperature.df, plot(Fold,-log10(FDR),
pch=20, main="Differential peak accessibility between 39 and 37",
col = ifelse(FDR<0.05 & abs(Fold)>1, "red", "black"),
cex = 0.75))
abline(v=c(-1,1), lty = "dotted")
abline(h=-log10(0.05), lty = "dotted")

## Adding transparency to the volcano plot
## It is hard to see the density of the points, so change transparency 
cols = c(rgb(255, 0, 0, max = 255, alpha = 25), rgb(0, 0, 0, max = 255, alpha = 25))
with(SN4741_temperature.df, plot(Fold,-log10(FDR),
pch=20, main="Differential peak accessibility between 39 and 37",
col = ifelse(FDR<0.05 & abs(Fold)>1, cols[1], cols[2]),
cex = 0.75))
abline(v=c(-1,1), lty = "dotted")
abline(h=-log10(0.05), lty = "dotted")

# How many differential peaks are there? Also view in IGV to confirm belief in peaks called DE 

## Add column to identify the peaks by the location, so can copy paste into IGV
SN4741_temperature.df$name <- paste0(SN4741_temperature.df$seqnames, ":", SN4741_temperature.df$start, "-", SN4741_temperature.df$end)

# Find peaks sig up and down in 39
up39 <- subset(SN4741_temperature.df, Fold > 1 & FDR <  0.05)
down39 <- subset(SN4741_temperature.df, Fold < -1 & FDR <  0.05)

dim(up39)
dim(down39)

# For viewing in IGV, sort by FC within the sig peaks to find the highest impacted regions 
sorted_up39 <- up39[order(-up39$Fold),]
sorted_down39 <- down39[order(-down39$Fold),]

# Write out and submit to GREAT to get GO terms of nearby genes
write.table(sorted_up39, "DAR_up39.txt", quote = F, row.names = F)
write.table(sorted_down39, "DAR_down39.txt", quote = F, row.names = F)

save.image("37v39.RData")
