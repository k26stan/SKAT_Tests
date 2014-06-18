##########################################################################################
## MrOSexome #############################################################################
##########################################################################################

# Tools
GATK_JAR=/projects/janssen/Tools/gatk2.7-2/GenomeAnalysisTK.jar
REF_FA=/projects/janssen/ref/ref.fa
VCF_TOOLS=/projects/janssen/Tools/vcftools_0.1.11/bin/vcftools
PLINK=/projects/janssen/Tools/plink-1.07-i686/plink 
GENE_TABLE=/home/kstandis/HandyStuff/GG-Gene_Names_DB.txt
GET_ALL_GENE_ANNOT_PY=/home/kstandis/RV/Get_All_Gene_Annot.py
GET_GENE_ANNOT_PY=/home/kstandis/RV/Get_Gene_Annot.py
SKAT_R=/home/kstandis/RV/GENE_SKAT.R
COMPILE_GENE_COORDS_R=/home/kstandis/RV/COMPILE_GENE_COORDS.R
GENE_MANH=/home/kstandis/RV/Gene_Manhattan.R
PLOT_SKAT_TESTS_R=/home/kstandis/RV/PLOT_SKAT_TESTS.R
GET_FUNC_VARS_R=/home/kstandis/RV/GET_FUNC_VARS.R

##########################################################################################
## MrOSexome #############################################################################
##########################################################################################

# Names/Paths
echo \### Defining Set Variables and Paths at `date` \###
DATE=20140506
HOME_DIR=/home/kstandis/MrOSexome

# Files
GENE_LIST=${HOME_DIR}/2012_GWAS_Unq.txt
IN_VCF=${HOME_DIR}/VAR_ONLY_VCFs/MERGE.SNP.vcf
ANNOTS=${HOME_DIR}/ANNOTATIONS/SNPS/SUMMARY_k26stan_190_novel_annotation.txt
PHENO_FILE=${HOME_DIR}/PHENOTYPE/Categ_Pheno.txt
COV_FILE=F

cd ${HOME_DIR}

DIRS=(${GENE_LIST//\// })
GENE_LIST_NAME=${DIRS[${#DIRS[@]} - 1]} # Set Path To Genotype List
DIRS=(${PHENO_FILE//\// })
PHENO=${DIRS[${#DIRS[@]} - 1]} # Get Name of Phenotype File

# Set up Directory for today's adventure
echo \### Moving to Base Directory at `date`: \###
ASSOC=${HOME_DIR}/${DATE}_${GENE_LIST_NAME%%.txt}_${PHENO%%.txt}
echo $ASSOC
#mkdir ${ASSOC}
#cp ${GENE_LIST} ${ASSOC}
NEW_GENE_LIST=${ASSOC}/Alt_${GENE_LIST_NAME}
GENE_COORDS=${ASSOC}/Gene_Coords.txt
#cd ${ASSOC}

##########################################################################################
## JnJ ###################################################################################
##########################################################################################

# Names/Paths
echo \### Defining Set Variables and Paths at `date` \###
DATE=20140429
HOME_DIR=/projects/janssen/RV_Rare_Variant

# Files
GENE_LIST=${HOME_DIR}/2013_Plenge_Vars.txt
IN_VCF=${HOME_DIR}/../VCFs/5-MERGE/HC_FULL.RG.SNP.vcf.gz
ANNOTS=${HOME_DIR}/../ANNOTATE/JnJ_121613_all_annotations.txt.gz
PHENO_FILE=${HOME_DIR}/PH-PHENOTYPES/B_CAT_52.txt
COV_FILE=${HOME_DIR}/PH-PHENOTYPES/COV_Bas.txt

cd ${HOME_DIR}

DIRS=(${GENE_LIST//\// })
GENE_LIST_NAME=${DIRS[${#DIRS[@]} - 1]} # Set Path To Genotype List
DIRS=(${PHENO_FILE//\// })
PHENO=${DIRS[${#DIRS[@]} - 1]} # Get Name of Phenotype File

# Set up Directory for today's adventure
echo \### Moving to Base Directory at `date`: \###
ASSOC=${HOME_DIR}/${DATE}_${GENE_LIST_NAME%%.txt}_${PHENO%%.txt}
echo $ASSOC
#mkdir ${ASSOC}
#cp ${GENE_LIST} ${ASSOC}
NEW_GENE_LIST=${ASSOC}/Alt_${GENE_LIST_NAME}
GENE_COORDS=${ASSOC}/Gene_Coords.txt
#cd ${ASSOC}

##########################################################################################
## Tristan's Kidney Data #################################################################
##########################################################################################

# Names/Paths
echo \### Defining Set Variables and Paths at `date` \###
DATE=20140506
HOME_DIR=/projects/ps-jcvi/projects/Kidney

# Files
GENE_LIST=${HOME_DIR}/Test_SKAT_Cap.txt
IN_VCF=${HOME_DIR}/hc-all.vcf
ANNOTS=${HOME_DIR}/SUMMARY_sga-upload_novel_annotation.txt
PHENO_FILE=${HOME_DIR}/Reject.txt
COV_FILE=F

cd ${HOME_DIR}

DIRS=(${GENE_LIST//\// })
GENE_LIST_NAME=${DIRS[${#DIRS[@]} - 1]} # Set Path To Genotype List
DIRS=(${PHENO_FILE//\// })
PHENO=${DIRS[${#DIRS[@]} - 1]} # Get Name of Phenotype File

# Set up Directory for today's adventure
echo \### Moving to Base Directory at `date`: \###
ASSOC=${HOME_DIR}/${DATE}_${GENE_LIST_NAME%%.txt}_${PHENO%%.txt}
echo $ASSOC
mkdir ${ASSOC}
cp ${GENE_LIST} ${ASSOC}
NEW_GENE_LIST=${ASSOC}/Alt_${GENE_LIST_NAME}
GENE_COORDS=${ASSOC}/Gene_Coords.txt
cd ${ASSOC}













##################################################################

for G in `cat ${NEW_GENE_LIST}`
do
	Rscript ${SKAT_R} ${ASSOC} ${G} ${PHENO_FILE} ${COV_FILE}
done

# Testin the New Python Script
python ${GET_ALL_GENE_ANNOT_PY} ${GENE_COORDS} ${ANNOTS} ${ASSOC}

# Manhattan Plots
Rscript ${GENE_MANH} ${NEW_GENE_LIST} ${COV_FILE}

# Plot SKAT Tests
Rscript ${PLOT_SKAT_TESTS_R} ${ASSOC} ${NEW_GENE_LIST}











###############################################################
#######################    END OF DOC    ######################
###############################################################