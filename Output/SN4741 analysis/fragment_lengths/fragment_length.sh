#!/bin/bash
#SBATCH --job-name=fragment_length
#SBATCH --time=2:00:00

module load samtools

# Use samtools and a horrible awk script to extract the mapping distance flag from the bams 
# Run from BAM folder

for file in *.bam;
	do
	samtools view $file | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > ${file%.bam}.fragment_length_count.txt
done

exit 0
