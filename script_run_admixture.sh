#!/bin/bash
#-------------------------------------------------------------------------------------------------------------
# MGG_GH5 - STRUCTURE ANALYSIS
# script_run_admixture.sh
#-------------------------------------------------------------------------------------------------------------
#SBATCH -J admixture
#SBATCH -o admixture_%j.o 
#SBATCH -e admixture_%j.e
#SBATCH -t 01:00:00 # execution time
#SBATCH --mem=1GB # memory needed
#SBATCH -c 4
#-------------------------------------------------------------------------------------------------------------
# Load modules
module load cesga/2020 admixture/1.3.0
#-------------------------------------------------------------------------------------------------------------
# Check number of arguments
if [ $# -ne 2 ]  
then
	echo "Usage: $0 <DATASET> <NPOP>" 
	exit 1
fi
# Set arguments
DATASET=$1 # Data name without file extension
NPOP=$2    # Number of populations
#-------------------------------------------------------------------------------------------------------------

### DIRECTORIES ###################

# Working directory
WDIR=$PWD

# Output directories
if [ ! -d results_admixture/$DATASET ] ; then mkdir -p results_admixture/$DATASET ; fi
cd results_admixture/$DATASET

### ADMIXTURE #####################

# Choosing the correct value for K (cross-validation error)

for (( K = 1; K <= $NPOP; K++ )); do
  admixture --cv $WDIR/data_pruned/$DATASET.bed $K | tee $DATASET_admixture_K${K}.out
done

# Outputs
# *.Q  The ancestry fractions 
# *.P  The allele frequencies of the inferred ancestral populations

# Extract the CV error for each corresponding K
awk '/CV/ {print $3,$4}' *.out | cut -c 4,7-20 > $DATASET.admixture.cv.error

