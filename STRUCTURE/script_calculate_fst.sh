#!/bin/bash
#-------------------------------------------------------------------------------------------------------------
# MGG_GH5 - STRUCTURE ANALYSIS
# script_calculate_fst.sh
#-------------------------------------------------------------------------------------------------------------
#SBATCH -J fst
#SBATCH -o fst_%j.o 
#SBATCH -e fst_%j.e
#SBATCH -t 01:00:00 # execution time
#SBATCH --mem=1GB # memory needed
#SBATCH -c 4
#-------------------------------------------------------------------------------------------------------------
# Load modules
module load cesga/2020 plink/1.9b5
#-------------------------------------------------------------------------------------------------------------
# Check number of arguments
if [ $# -ne 2 ]  
then
	echo "Usage: $0 <DATASET> <POPLIST>" 
	exit 1
fi
# Set arguments
DATASET=$1 # Data name without file extension
POPLIST=$2 # File with population assignment for each sample without file extension
#-------------------------------------------------------------------------------------------------------------

### DIRECTORIES ###################

# Working directory
WDIR=$PWD

# Output directories
if [ ! -d results_fst/$DATASET ] ; then mkdir -p results_fst/$DATASET ; fi
cd results_fst/$DATASET

### FST ###########################

plink --bfile $WDIR/data_pruned/$DATASET --fst --within $WDIR/$POPLIST.txt --out ${DATASET}_fst

# Outputs
# *.fst

