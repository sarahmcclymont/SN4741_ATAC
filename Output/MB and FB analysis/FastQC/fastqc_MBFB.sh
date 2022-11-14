#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --time=2:00:00
#SBATCH -N 1
#SBATCH -n 10

module load fastqc

## fastqc v0.11.9
fastqc -t 10 *.fastq.gz

# Generate report with MultiQC
## pip3 install --user multiqc #v1.13
multiqc .

# Clean up files
mkdir FastQC
mv *fastqc* FastQC

exit 0

```
