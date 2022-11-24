library(DiffBind) #v3.8.1

# Make the sample sheet containing all 37, 39, MB, and FB libraries
SampleID <- c("37-1", "37-2", "37-3", "37-4", "39-1", "39-2", "39-3", "39-4", "MB1", "MB2", "MB3", "FB1", "FB2", "FB3")
Condition <- c(rep("t37", 4), rep("t39", 4), rep("MB", 3), rep("FB", 3)) # add character in front of 37 and 39 because DiffBind doesn't love mask names that are numbers not characters in the future steps
Replicate <- c(rep(seq(1:4), 2), rep(seq(1:3), 2))
bamReads <- c("37-1_SN4741_S1.bam", "37-2_SN4741_S3.bam", "37-3_SN4741_S2.bam", "37-4_SN4741_S8.bam", "39-1_SN4741_S4.bam", "39-2_SN4741_S7.bam", "39-3_SN4741_S5.bam", "39-4_SN4741_S6.bam", "MB1_ATAC.bam", "MB2_ATAC.bam", "MB3_ATAC.bam", "FB1_ATAC.bam", "FB2_ATAC.bam", "FB3_ATAC.bam")
Peaks <- c("37-1_SN4741_S1_peaks_blacklisted.bed", "37-2_SN4741_S3_peaks_blacklisted.bed", "37-3_SN4741_S2_peaks_blacklisted.bed", "37-4_SN4741_S8_peaks_blacklisted.bed", "39-1_SN4741_S4_peaks_blacklisted.bed", "39-2_SN4741_S7_peaks_blacklisted.bed", "39-3_SN4741_S5_peaks_blacklisted.bed", "39-4_SN4741_S6_peaks_blacklisted.bed", "MB1_ATAC_peaks_blacklisted.bed", "MB2_ATAC_peaks_blacklisted.bed", "MB3_ATAC_peaks_blacklisted.bed", "FB1_ATAC_peaks_blacklisted.bed", "FB2_ATAC_peaks_blacklisted.bed", "FB3_ATAC_peaks_blacklisted.bed")
PeakCaller <- c(rep("bed", 14))
sampleData_comparison <- data.frame(cbind(SampleID, Condition, Replicate, bamReads, Peaks, PeakCaller))

## PCA and correlation heatmaps of different stages of data analysis
pdf("MB_comparison.pdf")

# Load the data into Diffbind
# Consensus peakset: peak must be present in ≥2 libraries
comparison <- dba(sampleSheet=sampleData_comparison, minOverlap = 2)

# No read counts - just peak presence/absence
dba.plotHeatmap(comparison) 
dba.plotPCA(comparison, DBA_CONDITION, label = DBA_ID)

# Count reads overlapping peaks and do not normalize
comparison_counts <- dba.count(comparison, summits = 100, bRemoveDuplicates = TRUE, score=DBA_SCORE_READS)
save.image() # count takes a long time, save.

# Plot heatmap and PCA considering raw read counts
dba.plotHeatmap(comparison_counts, score = DBA_SCORE_READS)
dba.plotPCA(comparison_counts, score = DBA_SCORE_READS, DBA_CONDITION, label = DBA_ID)

# Normalize the read counts
# Background, RLE normalization
comparison_normalized <- dba.normalize(comparison_counts, background = T, normalize = "RLE",  method=DBA_DESEQ2)

# Plot heatmap and PCA considering normalized read counts
dba.plotHeatmap(comparison_normalized, score = DBA_SCORE_NORMALIZED)
dba.plotPCA(comparison_normalized, score = DBA_SCORE_NORMALIZED, DBA_CONDITION, label = DBA_ID)


## Venn diagrams 

# Don't require a peak to be in multiple sets; it shows up in one library, we'll count it.
venn <- dba.peakset(comparison, consensus=DBA_CONDITION, minOverlap=1) 
dba.plotVenn(venn, venn$masks$Consensus)

## Upset plots 

library(UpSetR) #v1.4.0

# Just take the values from the Venn diagram - don't ask this package to find the overlaps 
# Will be annoying to type, read, and update in the future, but don't trust it to understand bed intervals. 
# Four way venn diagrams have 15 possible comparisons/compartments. Make sure you have them all. 

overlap_upset <- c(MB = 28957, 
                   FB = 14868, 
                   t37 = 24495, 
                   t39 = 3180, 
                   # Pairwise
                   `MB&FB` = 24479, 
                   `MB&t37` = 1455, 
                   `MB&t39` = 183,
                   `FB&t37` = 250,
                   `FB&t39` = 72,
                   `t37&t39` = 38182,
                   # Three way
                   `MB&FB&t37` = 2306,
                   `MB&FB&t39` = 297,
                   `MB&t37&t39` = 5060,
                   `FB&t37&t39` = 883,
                   # all 
                   `MB&FB&t37&t39` = 20667) 

upset(fromExpression(overlap_upset), order.by = "freq")
upset(fromExpression(overlap_upset))

dev.off()
