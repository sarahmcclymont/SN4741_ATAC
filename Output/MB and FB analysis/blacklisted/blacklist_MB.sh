#!/bin/bash
#SBATCH --job-name=blacklisted
#SBATCH --time=1:00:00

# Load modules
module load bedtools

# Define path to blacklist 
blacklist=/home/ext-smcclymont/SN4741_paper/mm10_ATAC_blacklist.bed

# Remove peak file header 
# Rearrange columns to bed format (chr start end name score strand), spoof strand column
  
# ! Check number of lines to remove from header and adjust code accordingly ! 
for file in *_peaks.xls;
	do
	awk 'NR > 26 { print }' < $file > temp
	awk -v OFS='\t' {'print $1,$2,$3,$10,$6,"+"'} temp > ${file%.xls}.bed
done

rm temp

# Count number of peaks called in each sample
wc -l *peaks.bed > peak_counts

# Remove the blacklisted regions from peaks and quantify number of peaks lost
# Use bedtools v2.30.0 to intersect peaks with blacklisted regions and exclude

for file in *peaks.bed; do
	bedtools intersect -v -wa -a $file -b $blacklist > ${file%.bed}_blacklisted.bed
done

wc -l *blacklisted* > blacklisted_counts

exit 0
