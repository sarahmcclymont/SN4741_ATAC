#!/bin/bash
#SBATCH --job-name=FRIP
#SBATCH --time=10:00:00
#SBATCH -N 1
#SBATCH -n 12

# Use deeptools v3.5.1 to plot enrichment of reads (bams) over the blacklisted peaks - once per lib
# Run from bams folder

module load libjpeg
module load bzip2

for file in *.bam; do
  plotEnrichment -b $file --BED ../blacklisted/${file%.bam}_peaks_blacklisted.bed -p 12 --outRawCounts ${file%.bam}.counts -v
done

exit 0

