#!/bin/bash
#-------------------------------------------------------------------------------------------------------------
# script_get_data.sh
#-------------------------------------------------------------------------------------------------------------
#SBATCH -J get_data
#SBATCH -o get_data_%j.o 
#SBATCH -e get_data_%j.e
#SBATCH -t 01:00:00 # execution time
#SBATCH --mem=16GB # memory needed
#SBATCH -c 4

module load cesga/2020 bcftools/1.19
#-------------------------------------------------------------------------------------------------------------

### WORKING DIRECTORY ###########################

WDIR=$PWD

### OUTPUT DIRECTORY ############################

if [ ! -d data ] ; then mkdir data ; fi
cd data

### 1000 GENOMES PROJECT ########################

# VARIABLES

N_POP=20 # number of individuals per population
POP_DIST=("YRI" "LWK" "GBR" "GIH" "CHB" "PEL" ) # case scenario 'geographically distant sampling' 
 # n_pop = 6
 # YRI: Yoruba, Africa
 # LWK: Luhya, Africa
 # GBR: British, Europe 
 # GIH: Gujarati, South Asia
 # CHB: Han Chinese, East Asia
 # PEL: Peruvian, America

POP_FINER=("GBR" "IBS" "TSI" "PJL" "BEB" "CDX" "CHB") # case scenario 'sampling with higher spatial resolution' 
 # n_pop = 7
 # GBR: British, Europe 
 # IBS: Iberian, Europe
 # TSI: Toscani, Europe
 # PJL: Punjabi, South Asia
 # BEB: Bengali, South Asia
 # CDX: Dai Chinese, East Asia
 # CHB: Han Chinese, East Asia

# DOWNLOAD COMPRESSED VCF AND PANEL WITH IDS AND POPULATION CODES

wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

# TABLE WITH POPULATION INFORMATION

awk 'NR>1{print $1"\t"$2"\t"$3}' integrated_call_samples_v3.20130502.ALL.panel > sample_pop.tsv

# FUNCTION TO PICK RANDOM INDIVIDUALS

pick_rnd_ids() {
  local pop=$1 # population id
  local n=$2 # number of individuals to pick
  awk -v p="$pop" '$2==p{print $1}' sample_pop.tsv | shuf -n "${n}"
}

# MAKE LIST PER CASE SCENARIO

for POP in "${POP_DIST[@]}"; do
  pick_rnd_ids "${POP}" "${N_POP}" | while read IID; do
    echo "${IID}" >> human_pop_dist.ids
    echo -e "0\t${IID}\t${POP}" >> poplist_human_pop_dist.txt
  done
done

for POP in "${POP_FINER[@]}"; do
  pick_rnd_ids "${POP}" "${N_POP}" | while read IID; do
    echo "${IID}" >> human_pop_finer.ids
    echo -e "0\t${IID}\t${POP}" >> poplist_human_pop_finer.txt
  done
done

# VCF SUBSET BY SAMPLES

bcftools view -S human_pop_dist.ids -Oz -o chr22_pop_dist.vcf.gz ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
bcftools index chr22_pop_dist.vcf.gz

bcftools view -S human_pop_finer.ids -Oz -o chr22_pop_finer.vcf.gz ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
bcftools index chr22_pop_finer.vcf.gz








