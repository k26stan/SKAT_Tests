## Parameters for Troubleshooting SKAT Pipeline ##
## June 12, 2014 ##
## Kristopher Standish ##


###########################################################
## 1 ## Manual Inputs (Osteo Data) ########################
###########################################################

# Names/Paths
echo \### Defining Set Variables and Paths at `date` \###
DATE=20140613b
HOME_DIR=/home/kstandis/MrOSexome
#HOME_DIR=/Users/kstandis/Data/MrOSexome

# Files
GENE_LIST=${HOME_DIR}/2012_GWAS_Tb2.txt
IN_VCF=${HOME_DIR}/VAR_ONLY_VCFs/MERGE.SNP.vcf
ANNOTS=${HOME_DIR}/ANNOTATIONS/SNPS/SUMMARY_k26stan_190_novel_annotation.txt
PHENO_FILE=${HOME_DIR}/PHENOTYPE/Quant_Pheno.txt
PHENO_TYPE=C
COV_FILE=F
COVS=`echo DAS_BL AGE SEX RF_ACPA BMI DRUG EE`
EIG_VEC=F
PC_COUNT=2
START_STEP=1

COVS=`echo "$COVS" | sed 's/ /QQQ/g'`

###########################################################
## REFER TO GENE_ASSOC.sh FOR THE REST ####################
###########################################################

## Run The Script
/home/kstandis/RV/Gene_Assoc.sh \
${DATE} \
${HOME_DIR} \
${GENE_LIST} \
${IN_VCF} \
${ANNOTS} \
${PHENO_FILE} \
${PHENO_TYPE} \
${COV_FILE} \
${COVS} \
${EIG_VEC} \
${PC_COUNT}




###########################################################
## 1 ## Manual Inputs (Janss Data) ########################
###########################################################


## Names/Paths
echo \### Defining Set Variables and Paths at `date` \###
DATE=20140616
HOME_DIR=/projects/janssen/RV_Rare_Variant

## Files
GENE_LIST=${HOME_DIR}/2014_Plenge.txt
IN_VCF=${HOME_DIR}/../VCFs/5-MERGE/HC_FULL.SNP.vcf
ANNOTS=${HOME_DIR}/../ANNOTATE/JnJ_121613_all_annotations.txt.gz
PHENO_FILE=${HOME_DIR}/../ASSOCIATION/PH-PHENOTYPES/LT8_DEL_28_BL.txt
PHENO_TYPE=C
COV_FILE=${HOME_DIR}/../ASSOCIATION/PH-PHENOTYPES/COV.txt
COVS=`echo DAS_BL AGE SEX`
EIG_VEC=${HOME_DIR}/../ASSOCIATION/EIGEN/HC_FULL.eigenvec
PC_COUNT=2
START_STEP=10

COVS=`echo "$COVS" | sed 's/ /QQQ/g'`

## Run The Script
/home/kstandis/RV/Gene_Assoc.sh \
${DATE} \
${HOME_DIR} \
${GENE_LIST} \
${IN_VCF} \
${ANNOTS} \
${PHENO_FILE} \
${PHENO_TYPE} \
${COV_FILE} \
${COVS} \
${EIG_VEC} \
${PC_COUNT} \
${START_STEP}
