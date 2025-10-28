#!/bin/bash
#-------------------------------------------------------------------------------------------------------------
# MGG_GH5 - STRUCTURE ANALYSIS
# script_data_filtering.sh
#-------------------------------------------------------------------------------------------------------------
#SBATCH -J datafilt
#SBATCH -o datafilt_%j.o 
#SBATCH -e datafilt_%j.e
#SBATCH -t 00:10:00 # execution time
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
DATASET=$1 # Data name without extension 'vcf.gz'
#-------------------------------------------------------------------------------------------------------------

### DIRECTORIES ###################

# Working directory
WDIR=$PWD

# Output directories
mkdir data_pruned

cd data_pruned

### DATA FILTERING ################

# Get plink files ('*.bed','*.fam','*.bim') keeping only biallelic variants (SNPs)
plink2 --vcf $WDIR/data/$DATASET.vcf.gz --make-bed --max-alleles 2 --snps-only --out $DATASET.SNPs

# Assign chromosome-and-position-based IDs to unnamed variants
plink2 --bfile $DATASET.SNPs --set-missing-var-ids @:# --make-bed --out $DATASET.SNPid
rm $DATASET.SNPs*

# Filter out variants/samples with more than 10 percent missing calls (--geno/--mind)
# Filter out variants with a minimum allele frequency < 0.05 (--maf)
plink2 --bfile $DATASET.SNPid --geno 0.1 --mind 0.1 --maf 0.05 --make-bed --out $DATASET.filtered
rm $DATASET.SNPid*

# Prune data for linkage disequilibrium
plink2 --bfile $DATASET.filtered --indep-pairwise 50 5 0.9 --out $DATASET.plink2
plink2 --bfile $DATASET.filtered --extract $DATASET.plink2.prune.in --make-bed --out ${DATASET}_pruned
rm $DATASET.filtered*
rm *.in
rm *.out
rm *.log


