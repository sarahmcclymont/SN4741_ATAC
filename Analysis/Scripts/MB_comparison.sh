#!/bin/bash

#SBATCH --job-name=DiffBind
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH -n 48

module load r
	# load R v4.2.2 -- for compatibility with a bunch of DiffBind dependencies
	
Rscript --save MB_comparison.R
  
exit 0

