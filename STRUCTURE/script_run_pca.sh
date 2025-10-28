#!/bin/bash
#-------------------------------------------------------------------------------------------------------------
# MGG_GH5 - STRUCTURE ANALYSIS
# script_run_pca.sh
#-------------------------------------------------------------------------------------------------------------
#SBATCH -J pca
#SBATCH -o pca_%j.o 
#SBATCH -e pca_%j.e
#SBATCH -t 01:00:00 # execution time
#SBATCH --mem=1GB # memory needed
#SBATCH -c 4
#-------------------------------------------------------------------------------------------------------------
# Load modules
module load cesga/2020 gcc/system plink/2.00a2.3
#-------------------------------------------------------------------------------------------------------------
# Check number of arguments
if [ $# -ne 1 ]  
then
	echo "Usage: $0 <DATASET>" 
	exit 1
fi
# Set arguments
DATASET=$1 # Data name without file extension
#-------------------------------------------------------------------------------------------------------------

### DIRECTORIES ###################

# Working directory
WDIR=$PWD

# Output directories
if [ ! -d results_pca/$DATASET ] ; then mkdir -p results_pca/$DATASET ; fi
cd results_pca/$DATASET

### PCA ###########################

plink2 --bfile $WDIR/data_pruned/$DATASET --pca 10 --out ${DATASET}_pca

# Outputs
# *.eigenval The eigenvalues from our analysis
# *.eigenvec The eigenvectors from our analysis


