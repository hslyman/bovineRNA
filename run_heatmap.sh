#!/bin/bash
#SBATCH --time=30
#SBATCH --mem=8000

module load R

Rscript ./bovine_RNASeq_heatmap.R


