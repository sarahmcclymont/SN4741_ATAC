#!/bin/bash
#SBATCH --job-name=alignment
#SBATCH --time=30:00:00
#SBATCH -N 1
#SBATCH -n 12

# Load modules
  module load bowtie2 #v2.4.1
  module load samtools #v1.15.1 
  module load picard #v2.26.11

# Define path to genome indices 
  genome_index=/data/amccall2/work-zfs/NBB/SN4741-pcHiC/bowtie2_index/mm10

#Hard to make for loop when some samples have been sequenced twice but most just once. Exclude them from the loop by temporarily moving to a new dir; analyse separately
  
for file in *R1.fastq.gz; do

# Make life easier by making the file names simpler and into an extensionlesss variable 
  sample=${file%_R1.fastq.gz} 
  
# Write out sample name so can track progress/errors in slurm.out file 
  echo $sample

# Align with Bowtie2
  # Local alignment to softclip the ends
  # Increase pair distance to 1000 to avoid the 500bp blip on fragment length plot
  # Write out alignment stats (stderr) to file for multiqc summary
  bowtie2 -p 12 --local -X 1000 -x $genome_index -1 $file -2 ${sample}_R2.fastq.gz -S ${sample}.sam 2>${sample}.align.stats
  
# Convert to bam and sort with samtools
  # view -bu turns into uncompressed bam for easy piping into sort
  samtools view -bu ${sample}.sam | samtools sort -o ${sample}_sorted.bam
  
# Deduplicate with picard tools
  picard MarkDuplicates I=${sample}_sorted.bam O=${sample}.bam M=${sample}_markdup.stats 
  
# Index with samtools 
  samtools index ${sample}.bam
  
# Remove intermediate files 
#	rm ${sample}.sam
#	rm ${sample}_sorted.bam
  
# Call peaks!
	macs3 callpeak --nomodel -B -f BAMPE -t ${sample}.bam -n $sample --nolambda --gsize mm 
	
	done
	
exit 0
