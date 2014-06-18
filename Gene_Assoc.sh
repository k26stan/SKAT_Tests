## Notes about GENE_ASSOC.sh

###########################################################
## 1 ## Set up Paths ######################################
###########################################################
echo \### 1 - `date` \###
echo \### Defining Set Variables and Paths \###

###########################################################
## Manually Input Parameters ##

# Names/Paths
DATE=$1
HOME_DIR=$2

# Parameters and Files
GENE_LIST=$3 # List of Genes you want to test
IN_VCF=$4 # Path to VCF File with all the Variants
ANNOTS=$5 # Path to Annotation File
PHENO_FILE=$6 # Which Phenotype File are you using?
PHENO_TYPE=$7 # Is phenotype (B)inary or (C)ontinuous?
COV_FILE=$8 # Path to Covariate File or "F"
COVS=$9 # Which Covariates to Include?
EIG_VEC=${10} # Output from Plink's --pca command (MAF>1%) of "F"
PC_COUNT=${11} # How many PCs to Include as Covariates?
START_STEP=${12}

###########################################################
## Constant Paths ##

# Tools
GATK_JAR=/projects/janssen/Tools/gatk2.7-2/GenomeAnalysisTK.jar
REF_FA=/projects/janssen/ref/ref.fa
VCF_TOOLS=/projects/janssen/Tools/vcftools_0.1.11/bin/vcftools
PLINK=/projects/janssen/Tools/plink_linux_x86_64/plink 
GENE_TABLE=/home/kstandis/HandyStuff/GG-Gene_Names_DB.txt
GET_ALL_GENE_ANNOT_PY=/home/kstandis/RV/Get_All_Gene_Annot.py
GET_GENE_ANNOT_PY=/home/kstandis/RV/Get_Gene_Annot.py
SKAT_R=/home/kstandis/RV/GENE_SKAT.R
COMPILE_GENE_COORDS_R=/home/kstandis/RV/COMPILE_GENE_COORDS.R
GENE_MANH=/home/kstandis/RV/Gene_Manhattan.R
PLOT_SKAT_TESTS_R=/home/kstandis/RV/PLOT_SKAT_TESTS.R
GET_FUNC_VARS_R=/home/kstandis/RV/GET_FUNC_VARS.R
MAKE_COV_TAB_R=/home/kstandis/RV/Make_Cov_Tab.R

###########################################################
## Pull some Info out of Parameters ##

# Get Names of Specific Files
DIRS=(${GENE_LIST//\// })
GENE_LIST_NAME=${DIRS[${#DIRS[@]} - 1]} # Get Name of Genotype List
DIRS=(${PHENO_FILE//\// })
PHENO=${DIRS[${#DIRS[@]} - 1]} # Get Name of Phenotype File

# Specify list of Covariates to include (for command ad for filename)
if [ $PC_COUNT -eq 0 ] ; then
COVS_COMMAND=`echo "${COVS}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${COVS}" | sed 's/QQQ/_/g'`
else
PCS=`seq 1 ${PC_COUNT}`
PCS_COMMAND=`echo "PC"${PCS} | sed 's/ /QQQPC/g'`
COVS_COMMAND=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/_/g'`
fi

if [ $COV_FILE = "F" ] ; then
COVS_COMMAND=`echo "${PCS_COMMAND}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${PCS_COMMAND}" | sed 's/QQQ/_/g'`
fi

# Specify commands and extensions for Cont vs Bin Phenotype
if [ $PHENO_TYPE = "C" ] ; then
SUFFIX=linear
else
SUFFIX=logistic
fi

# Set up Directory for today's adventure
echo \### Moving to Base Directory at `date`: \###
echo $ASSOC
ASSOC=${HOME_DIR}/${DATE}_${GENE_LIST_NAME%%.txt}_${PHENO%%.txt}_${COVS_FILENAME}
mkdir ${ASSOC}
cp ${GENE_LIST} ${ASSOC}
NEW_GENE_LIST=${ASSOC}/Alt_${GENE_LIST_NAME}
GENE_COORDS=${ASSOC}/Gene_Coords.txt
NEW_COV_FILE=${ASSOC}/Cov_w_PCs.txt
cd ${ASSOC}

# Specify a File to which to Write Updates
UPDATE_FILE=${ASSOC}/Update.txt

if [ "$START_STEP" -le 1 ]; then
echo `date` "1 - Got all the paths taken care of...I think" > ${UPDATE_FILE}
echo "This should be the GATK directory:" ${GATK_JAR} >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 2 ## Make sure I've got coordinates ####################
###########################################################
if [ "$START_STEP" -le 2 ]; then
 # Either give coordinates in input or search for them
   # Accommodate either situation
   # Giving coordinates will end up being the more reliable option
echo \### 2 - `date` \###
echo \### Retreiving and Compiling Gene Coordinates \###
echo `date` "2 - Retreiving and Compiling Gene Coordinates" >> ${UPDATE_FILE}

N_COL=`head -1 ${GENE_LIST} | awk '{print NF}'`
if [ $N_COL = 1 ] # If Only Gene Name is Given (Not Coordinates)
then
echo Only Gene Name Given
for G in `cat ${GENE_LIST}`
do
mkdir ${ASSOC}/${G}
awk -F $'\t' -v GENE="$G" '{ if ( $15 == GENE ) print $2":"$3"-"$4 }' ${GENE_TABLE} > ${ASSOC}/${G}/${G}.list
sed -i 's/chr//g' ${ASSOC}/${G}/${G}.list
done
else # If Coordinates are Given Too
echo Gene Names and Coordinates Provided
IFSo=$IFS
IFS=$'\n' # Makes it so each line is read whole (not separated by tabs)
for LINE in `cat ${GENE_LIST}` # Make Coordinate Lists from Provided Info
do
G="$(echo "$LINE" | awk '{print $1}')" # Gene Name is 1st Tab-Delimited Column
COORDS=`echo "${LINE}" | awk '{print $2}'` # Coordinates are 2nd Tab-Delimited Column (CHR:START-STOP)
mkdir ${ASSOC}/${G} # Make Directory from Gene Name
echo ${COORDS} > ${ASSOC}/${G}/${G}.list
done
cat ${GENE_LIST} | awk '{print $1}' > ${ASSOC}/${GENE_LIST_NAME}
GENE_LIST=${ASSOC}/${GENE_LIST_NAME}
IFS=$IFSo # Reset
fi

echo `date` "Pulled out all the Gene Coordinates" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 3 ## Compile Gene Coordinates ##########################
###########################################################
if [ "$START_STEP" -le 3 ]; then
 # Rscript - Compile_Gene_Coords.R ${ASSOC} ${GENE_LIST}
   # Outputs "Gene_Coords.txt" file
   # Outputs $NEW_GENE_LIST = ".../Alt_TEST_GENES.txt" file (which omits Genes with missing coordinates)
   # Re-writes all ${G}/${G}.list files with an extra 10K on either end
echo \### 3 - `date` \###
echo \### Compiling Gene Coordinates and adding 10K to each end - Rscript \###
echo `date` "3 - Compiling Gene Coordinates and adding 10K to each end - Rscript" >> ${UPDATE_FILE}

Rscript ${COMPILE_GENE_COORDS_R} ${ASSOC} ${GENE_LIST}

echo `date` "Compiled Gene Coordinates and added 10K" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 4 ## Pull out SNPs for each Gene #######################
###########################################################
if [ "$START_STEP" -le 4 ]; then
 # Use GATK b/c it'll end up being faster
   # Save each to it's own directory/file
echo \### 4 - `date` \###
echo \### Pull of Variants for each Gene - GATK \###
echo `date` "4 - Pull of Variants for each Gene - GATK" >> ${UPDATE_FILE}

for G in `cat ${NEW_GENE_LIST}`
do
	java -Xmx3g -jar ${GATK_JAR} \
	-T SelectVariants \
	-R ${REF_FA} \
	-V ${IN_VCF} \
	-L ${ASSOC}/${G}/${G}.list \
	-o ${ASSOC}/${G}/${G}.vcf
	echo \# Finished ${G} at `date` \#
	echo `date` "GATK -" ${G} >> ${UPDATE_FILE}
done

echo `date` "Done pulling out variants using GATK" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 5 ## Get Gene Annotations from Cypher/SG-ADVISER File ##
###########################################################
if [ "$START_STEP" -le 5 ]; then
 # Python - Get_All_Gene_Annot.py
echo \### 5 - `date` \###
echo \### Pull of Annotations for each Gene - Python \###
echo `date` "5 - Pull of Annotations for each Gene - Python" >> ${UPDATE_FILE}

python ${GET_ALL_GENE_ANNOT_PY} ${GENE_COORDS} ${ANNOTS} ${ASSOC} # Output is Annotation File for each

echo `date` "Done pulling out Annotations using Python" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 6 ## Pull out Functional Variants from Annotations #####
###########################################################
if [ "$START_STEP" -le 6 ]; then
 # Rscript - Get_Func_Vars.R ${ASSOC} ${NEW_GENE_LIST} ${ANNOTS}
   # Pull out functional variants from Annotation File
echo \### 6 - `date` \###
echo \### Identify Functional Variants - Rscript \###
echo `date` "6 - Identify Functional Variants - Rscript" >> ${UPDATE_FILE}

Rscript ${GET_FUNC_VARS_R} ${ASSOC} ${NEW_GENE_LIST} ${ANNOTS}

echo `date` "Done Identifying Functional Variants using Rscript" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 7 ## Re-Genotype Missing Variants to 0/0 ###############
###########################################################
if [ "$START_STEP" -le 7 ]; then
echo \### 7 - `date` \###
echo \### Regenotype Missing SNPs as 0/0 \###
echo `date` "7 - Regenotype Missing SNPs as 0/0" >> ${UPDATE_FILE}

# Regenotype VCF files from ./. to 0/0
for G in `cat ${NEW_GENE_LIST}`
do
	sed -i 's/\.\/\./0\/0/g' ${ASSOC}/${G}/${G}.vcf
done

echo `date` "Regenotyping Done" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 8 ## Convert files to Plink format #####################
###########################################################
if [ "$START_STEP" -le 8 ]; then
 # VCFtools - "Plink" command
echo \### 8 - `date` \###
echo \### Converting Files to PED/01 Format - VCFtools \###
echo `date` "8 - Converting Files to PED/01 Format - VCFtools" >> ${UPDATE_FILE}

for G in `cat ${NEW_GENE_LIST}`
do
	mkdir ${ASSOC}/${G}/PLINK_FILES
	$VCF_TOOLS --vcf ${ASSOC}/${G}/${G}.vcf --plink --out ${ASSOC}/${G}/PLINK_FILES/${G}
	${VCF_TOOLS} --vcf ${ASSOC}/${G}/${G}.vcf --012 --out ${ASSOC}/${G}/${G}
done

echo `date` "Converting VCFs to PED/01 Format Done" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 9 ## Make New Covariate File Containing PCs ############
###########################################################
if [ "$START_STEP" -le 9 ]; then
 # If no PCs exist, make them
   # Convert to PED and then calculate PCs
echo \### 9 - `date` \###
echo \### Compile Covariates with PCs \###
echo `date` "9 - Compile Covariates with PCs" >> ${UPDATE_FILE}

# If No Principal Components Exist, Make Them
if [ $EIG_VEC = "F" ] ; then
# Convert Original VCF to PED File
${VCF_TOOLS} --vcf ${IN_VCF} --plink --out ${ASSOC}/Full
# Use PED to run PCA
${PLINK} --file ${ASSOC}/Full \
--pca header \
--allow-no-sex \
--out ${ASSOC}/Full
EIG_VEC=${ASSOC}/Full.eigenvec
fi

# Make Covariate File for this Run (That includes PCs)
if [ $COV_FILE = "F" ] ; then
cp ${EIG_VEC} ${NEW_COV_FILE}
else
Rscript ${MAKE_COV_TAB_R} ${EIG_VEC} ${COV_FILE} ${NEW_COV_FILE}
fi

echo `date` "Compiling Covariates Done" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 10 ## Run Single-Locus Association on each Variant #####
###########################################################
if [ "$START_STEP" -le 10 ]; then
 # Plink - Logistic/Linear
echo \### 10 - `date` \###
echo \### Running Single-Locus Association - PLINK \###
echo `date` "10 - Running Single-Locus Association - PLINK" >> ${UPDATE_FILE}

for G in `cat ${NEW_GENE_LIST}`
do
${PLINK} --file ${ASSOC}/${G}/PLINK_FILES/${G} \
--pheno ${PHENO_FILE} \
--covar ${NEW_COV_FILE} --covar-name ${COVS_COMMAND} \
--${SUFFIX} hide-covar --adjust \
--allow-no-sex \
--maf 0.01 \
--out ${ASSOC}/${G}/PLINK_FILES/PLINK_${G}
done

echo `date` "Single-Locus Associations Done" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 11 ## Make Manhattan Plot for Single Locus Tests #######
###########################################################
if [ "$START_STEP" -le 11 ]; then
 # Rscript
 # Also compiles list of SNPs beyond a given threshold (0.1?)
echo \### 11 - `date` \###
echo \### Make Manhattan Plot of Genes - Rscript \###
echo `date` "11 - Make Manhattan Plot of Genes - Rscript" >> ${UPDATE_FILE}

Rscript ${GENE_MANH} ${NEW_GENE_LIST} ${COV_FILE} ${PHENO_TYPE}

echo `date` "Plotting the Skyline Done" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 12 ## Run SKAT Tests ###################################
###########################################################
if [ "$START_STEP" -le 12 ]; then
 # Rscript
 # Use different collapsed sets
 # Different Tests
echo \### 12 - `date` \###
echo \### Run SKAT Tests - Rscript \###
echo `date` "12 - Run SKAT Tests - Rscript" >> ${UPDATE_FILE}

Rscript ${SKAT_R} ${ASSOC} ${NEW_GENE_LIST} ${PHENO_FILE} ${PHENO_TYPE} ${NEW_COV_FILE} ${COVS_COMMAND}

echo `date` "SKAT Tests Done" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 13 ## Plot Results of SKAT Tests #######################
###########################################################
if [ "$START_STEP" -le 13 ]; then
 # Rscript
 # TBD...
echo \### 13 - `date` \###
echo \### Plot SKAT Results - Rscript \###
echo `date` "13 - Plot SKAT Results - Rscript" >> ${UPDATE_FILE}

Rscript ${PLOT_SKAT_TESTS_R} ${ASSOC}

echo `date` "Plotting SKAT Results Done" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 14 ## END OF DOC...SO FAR ##############################
###########################################################
echo `date` "DONE!!!" >> ${UPDATE_FILE}
echo \############################################ >> ${UPDATE_FILE}
printf "DONE!!!"









