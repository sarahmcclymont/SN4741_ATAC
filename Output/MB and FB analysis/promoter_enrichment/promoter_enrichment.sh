#!/bin/bash
#SBATCH --job-name=MB_enrich_proms
#SBATCH --time=15:00:00
#SBATCH -N 1
#SBATCH -n 12

# Use deeptools v3.5.1 to plot enrichment of reads (bams) over the promoters
# Run from bams folder

module load libjpeg
module load bzip2

promoters=/home/ext-smcclymont/SN4741_paper/mm10_promoters_1000bp.bed

for file in *.bam; do
  sample=${file:0:3} # Just keeps the sample name with none of the extra ("SN4741_S#) that has been carrying through; keep it clean
	bamCoverage -p 12 -b $file -o ${sample}.bw
	computeMatrix reference-point -p 12 --referencePoint center -b 1000 -a 1000 -bs 50 -S ${sample}.bw -R $promoters -out ${sample}.matrix
	plotHeatmap --matrixFile ${sample}.matrix --outFileName ${sample}.heatmap.pdf --heatmapHeight 5 --heatmapWidth 5 
done

exit 0
