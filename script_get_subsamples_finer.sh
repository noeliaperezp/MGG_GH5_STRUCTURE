#!/bin/bash
#-------------------------------------------------------------------------------------------------------------
# MGG_GH5 - STRUCTURE ANALYSIS
# script_get_subsamples_finer.sh
#-------------------------------------------------------------------------------------------------------------
#SBATCH -J get_data
#SBATCH -o get_data_%j.o 
#SBATCH -e get_data_%j.e
#SBATCH -t 00:10:00 # execution time
#SBATCH --mem=16GB # memory needed
#SBATCH -c 4

module load cesga/2020 bcftools/1.19
#-------------------------------------------------------------------------------------------------------------

### WORKING DIRECTORY ###########################

WDIR=$PWD


### OUTPUT DIRECTORY ############################

if [ ! -d data_subsamples ] ; then mkdir data_subsamples ; fi
cd data_subsamples/


### POPULATIONS #################################

N_POP_FINER=50                                        # number of individuals per population
POP_FINER=("IBS" "TSI" "GBR" "FIN" "GIH")             # case scenario 'sampling with higher spatial resolution' 


### TABLE WITH POPULATION INFORMATION ###########

awk 'NR>1{print $1"\t"$2"\t"$3}' ../data/integrated_call_samples_v3.20130502.ALL.panel > sample_pop.tsv


### FUNCTION TO PICK RANDOM INDIVIDUALS #########

pick_rnd_ids() {
  local pop=$1 # population id
  local n=$2 # number of individuals to pick
  awk -v p="$pop" '$2==p{print $1}' sample_pop.tsv | shuf -n "${n}"
}


### MAKE LIST PER CASE SCENARIO #################

rm human_pop_finer.ids
rm poplist_human_pop_finer.txt

for POP in "${POP_FINER[@]}"; do
  pick_rnd_ids "${POP}" "${N_POP_FINER}" | while read IID; do
    echo "${IID}" >> human_pop_finer.ids
    echo -e "0\t${IID}\t${POP}" >> poplist_human_pop_finer.txt
  done
done


### VCF SUBSET BY SAMPLES #######################

bcftools view -S human_pop_finer.ids -Oz -o chr22_pop_finer.vcf.gz ../data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
bcftools index chr22_pop_finer.vcf.gz


