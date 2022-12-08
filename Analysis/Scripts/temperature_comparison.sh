#!/bin/bash

#SBATCH --job-name=temperature_comp
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH -n 48

module load r

Rscript --save temperature_comparison.R
  # loads R v4.2.2 -- for compatibility with a bunch of DiffBind dependencies

exit 0

