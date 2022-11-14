#!/bin/bash
#SBATCH --job-name=reseq_alignment
#SBATCH --time=30:00:00
#SBATCH -N 1
#SBATCH -n 12

# Load modules
  module load bowtie2 #v2.4.1
  module load samtools #v1.15.1 
  module load picard #v2.26.11

# Define path to genome indices 
  genome_index=/data/amccall2/work-zfs/NBB/SN4741-pcHiC/bowtie2_index/mm10
  
# Align individually
  bowtie2 -p 12 --local -X 1000 -x $genome_index -1 FB2_ATAC_R1.fastq.gz,FB2_reseq_ATAC_R1.fastq.gz -2 FB2_ATAC_R2.fastq.gz,FB2_reseq_ATAC_R2.fastq.gz -S FB2_ATAC.sam 2>FB2.align.stats
  bowtie2 -p 12 --local -X 1000 -x $genome_index -1 MB3_ATAC_R1.fastq.gz,MB3_reseq_ATAC_R1.fastq.gz -2 MB3_ATAC_R2.fastq.gz,MB3_reseq_ATAC_R2.fastq.gz -S MB3_ATAC.sam 2>MB3.align.stats

# Reinitiate loop
for file in *.sam; do

# Make life easier by making the file names simpler and into an extensionlesss variable 
  sample=${file%.sam} 
  
# Write out sample name so can track progress/errors in slurm.out file 
  echo $sample
  
# Convert to bam and sort with samtools
  # view -bu turns into uncompressed bam for easy piping into sort
  samtools view -bu $file | samtools sort -o ${sample}_sorted.bam
  
# Deduplicate with picard tools
  picard MarkDuplicates I=${sample}_sorted.bam O=${sample}.bam M=${sample}_markdup.stats 
  
# Index with samtools 
  samtools index ${sample}.bam
  
# Remove intermediate files 
#	rm ${sample}.sam
#	rm ${sample}_sorted.bam
  
# Call peaks!
	macs3 callpeak --nomodel -B -f BAMPE -t ${sample}.bam -n $sample --nolambda --gsize mm 

# Remove intermediate files 
#	rm ${sample}.sam
#	rm ${sample}_sorted.bam
	
	done
	
exit 0
