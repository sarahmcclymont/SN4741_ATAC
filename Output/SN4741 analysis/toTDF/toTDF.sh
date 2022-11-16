#!/bin/bash
#SBATCH --job-name=toTDF
#SBATCH --time=10:00:00

module load Java

for file in *.bdg;
	do
	~/Software/IGV_2.15.2/igvtools toTDF $file ${file%.bdg} mm10
done

exit 0
