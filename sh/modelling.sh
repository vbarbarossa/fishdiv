#!/bin/bash
#SBATCH -N 1
#SBATCH -t 2-00:00:00
#SBATCH -p normal
#SBATCH --output=sh/modelling.out
#SBATCH --mail-type=END
#SBATCH --mail-user=v.barbarossa@cml.leidenuniv.nl

module load 2019
module load R/3.5.1-foss-2018b

Rscript R/3_modelling_untransformed.R
