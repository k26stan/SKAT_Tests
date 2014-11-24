## Gene_Assoc.sh
## Shell Script to Perform Rare Variant SKAT Tests on Specific Genes ##
## This takes in specific parameters and calls several different modules to complete the tests ##
## Originally Developed April/June, 2014 ##
## Refined August, 2014 ##
## Kristopher Standish ##
## UCSD/JCVI ##

###########################################################
## 1 ## Set up Paths/Parameters/etc... ####################
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
# IN_VCF=$4 # Path to VCF File with all the Variants
VAR_FILE=$4 # Maybe switch this from VCF to whatever somebody wants
ANNOTS=$5 # Path to Annotation File
PHENO_FILE=$6 # Which Phenotype File are you using?
PHENO_TYPE=$7 # Is phenotype (B)inary or (C)ontinuous?
COV_FILE=$8 # Path to Covariate File or "F"
COVS=$9 # Which Covariates to Include?
EIG_VEC=${10} # Output from Plink's --pca command (MAF>1%) of "F"
PC_COUNT=${11} # How many PCs to Include as Covariates?
START_STEP=${12} # Which Step do you want to start on?
GENE_START=${13} # When running SKAT tests, which Gene to start on (Integer: 1-??)

###########################################################
## Constant Paths ##
RV_DIR=/home/kstandis/RV

# Public Tools
GATK_JAR=${RV_DIR}/Tools/GenomeAnalysisTK.jar
REF_FA=${RV_DIR}/Tools/ref.fa
VCF_TOOLS=${RV_DIR}/Tools/vcftools
PLINK=${RV_DIR}/Tools/plink 
GENE_TABLE=${RV_DIR}/Tools/GG-Gene_Names_DB.txt
# Custom Scripts
COMPILE_GENE_COORDS_R=${RV_DIR}/Compile_Gene_Coords.R
GET_ALL_GENE_ANNOT_PY=${RV_DIR}/Get_All_Gene_Annot.py
GET_FUNC_VARS_R=${RV_DIR}/Get_Func_Vars.R
MAKE_COV_TAB_R=${RV_DIR}/Make_Cov_Tab.R
GENE_MANH=${RV_DIR}/Gene_Manhattan.R
SKAT_R=${RV_DIR}/Gene_SKAT.R
PLOT_SKAT_TESTS_R=${RV_DIR}/Plot_SKAT_Tests.R

###########################################################
## Pull some Info out of Parameters ##

# Get Names of Specific Files
DIRS=(${GENE_LIST//\// })
GENE_LIST_NAME=${DIRS[${#DIRS[@]} - 1]} # Get Name of Geno List
DIRS=(${PHENO_FILE//\// })
PHENO=${DIRS[${#DIRS[@]} - 1]} # Get Name of Phenotype File

echo $PHENO
echo $DIRS
echo $GENE_LIST_NAME
echo whatever

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

# Incorporate Country/Site of Study as Binary Covariate (if Included)
if [[ $COVS == *COUN* ]]
then
COVS_COMMAND=`echo $COVS_COMMAND | sed 's/COUN/CN_ARG,CN_AUS,CN_COL,CN_HUN,CN_LTU,CN_MEX,CN_MYS,CN_NZL,CN_POL,CN_RUS,CN_UKR/g'`
fi
echo Covariates \(Command\): ${COVS_COMMAND}

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
if [ ! -d "${ASSOC}" ] ; then mkdir ${ASSOC} ; fi
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
## 2 ## Verify/Retrieve Gene Coordinates ##################
###########################################################
if [ "$START_STEP" -le 2 ]; then
 # Either give coordinates in input or search for them
   # Accommodate either situation
   # Giving coordinates will end up being the more reliable option
echo \### 2 - `date` \###
echo \### Retreiving and Compiling Gene Coordinates \###
echo `date` "2 - Retreiving and Compiling Gene Coordinates" >> ${UPDATE_FILE}

# Determine How Many Columns are in Gene List
N_COL=`head -1 ${GENE_LIST} | awk '{print NF}'`

# If Only Gene Name is Given (Not Coordinates)
if [ $N_COL = 1 ]
then
echo Only Gene Name Given
for G in `cat ${GENE_LIST}`
do
if [[ $G == */* ]]
then
G=`echo ${G} | sed 's@/@_@g'`
fi
if [ ! -d "${ASSOC}/${G}" ] ; then mkdir ${ASSOC}/${G} ; fi
awk -F $'\t' -v GENE="$G" '{ if ( $15 == GENE ) print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$15 }' ${GENE_TABLE} > ${ASSOC}/${G}/${G}.exon.list
# awk -F $'\t' -v GENE="$G" '{ if ( $15 == GENE ) print $2":"$3"-"$4 }' ${GENE_TABLE} > ${ASSOC}/${G}/${G}.list
cat ${ASSOC}/${G}/${G}.exon.list | awk '{print $1":"$2"-"$3}' > ${ASSOC}/${G}/${G}.1.list
sed -i 's/chr//g' ${ASSOC}/${G}/${G}.1.list
done
# Else (if Coordinates are Given Also)
else
echo Gene Names and Coordinates Provided
IFSo=$IFS
IFS=$'\n' # Makes it so each line is read whole (not separated by tabs)
for LINE in `cat ${GENE_LIST}` # Make Coordinate Lists from Provided Info
do
G="$(echo "$LINE" | awk '{print $1}')" # Gene Name is 1st Tab-Delimited Column
COORDS=`echo "${LINE}" | awk '{print $2}'` # Coordinates are 2nd Tab-Delimited Column (CHR:START-STOP)
if [ ! -d "${ASSOC}/${G}" ] ; then mkdir ${ASSOC}/${G} ; fi # Make Directory from Gene Name
echo ${COORDS} > ${ASSOC}/${G}/${G}.1.list
awk -F $'\t' -v GENE="$G" '{ if ( $15 == GENE ) print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$15 }' ${GENE_TABLE} > ${ASSOC}/${G}/${G}.exon.list
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
   # Re-writes all ${G}/${G}.list files with an extra 5K on either end
echo \### 3 - `date` \###
echo \### Compiling Gene Coordinates and adding 5K to each end - Rscript \###
echo `date` "3 - Compiling Gene Coordinates and adding 5K to each end - Rscript" >> ${UPDATE_FILE}

Rscript ${COMPILE_GENE_COORDS_R} ${ASSOC} ${GENE_LIST}

echo `date` "Compiled Gene Coordinates and added 5K" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 4 ## Pull out SNPs for each Gene #######################
###########################################################
if [ "$START_STEP" -le 4 ]; then
 # Use GATK b/c it'll end up being faster
   # Save each to it's own directory/file
echo \### 4 - `date` \###
echo \### Pull Variants for each Gene \###
echo `date` "4 - Pull Variants for each Gene" >> ${UPDATE_FILE}

###########################################################
## Pull out variants for all genes of interest
 # To make a more manageable variant file for quicker access
echo \### 4a - `date` \###
echo \### Pull Variants for Genes into 1 File \###
echo `date` "4a - Pull Variants for Genes into 1 File" >> ${UPDATE_FILE}

## Pull out Genes into a Single VCF file w/ all Gene Variants
if [ ${VAR_FILE: -4} == ".vcf" -o ${VAR_FILE: -7} == ".vcf.gz" ] ; then
	# Pull out Variants in each Gene
	if [ ${VAR_FILE: -4} == ".vcf" ] ; then
		${VCF_TOOLS} \
		--vcf ${VAR_FILE} \
		--extract range ${ASSOC}/Gene_Coords.pl.list \
		--recode \
		--out ${ASSOC}/Var_File_1
	fi
	if [ ${VAR_FILE: -7} == ".vcf.gz" ] ; then
		${VCF_TOOLS} \
		--gzvcf ${VAR_FILE} \
		--extract range ${ASSOC}/Gene_Coords.pl.list \
		--recode \
		--out ${ASSOC}/Var_File_1
	fi
	# Re-Genotype VCF file for each gene
	sed -i 's/\.\/\./0\/0/g' ${ASSOC}/${G}/${G}.vcf
	# Convert to .ped file
	${VCF_TOOLS} \
	--vcf ${ASSOC}/Var_File_1 \
	--plink \
	--out ${ASSOC}/Var_File_2
	# Convert to .bed file
	${PLINK} \
	--vcf ${ASSOC}/Var_File_2 \
	--silent \
	--make-bed \
	--out ${ASSOC}/Var_File_3
	echo \# Pulled Vars for All Genes to 1 file - `date` \#
	echo `date` "All Genes Pulled to 1 File" >> ${UPDATE_FILE}
fi

## If File is ".ped" or ".bed", use Plink to pull out variants for all gene
 # Output should be .bed file for each gene
if [ ${VAR_FILE: -4} == ".ped" ] ; then
	# Output single variant BED file for All genes
	${PLINK} \
	--file ${VAR_FILE%%.ped} \
	--silent \
	--extract range ${ASSOC}/Gene_Coords.pl.list \
	--hardy midp \
	--make-bed \
	--memory 8000 \
	--threads 1 \
	--out ${ASSOC}/Var_File_3
	echo \# Pulled Vars for All Genes to 1 file -  `date` \#
	echo `date` "All Genes Pulled to 1 File" >> ${UPDATE_FILE}
fi
if [ ${VAR_FILE: -4} == ".bed" ] ; then
	# Output single variant BED file for All genes
	${PLINK} \
	--bfile ${VAR_FILE%%.bed} \
	--silent \
	--extract range ${ASSOC}/Gene_Coords.pl.list \
	--hardy midp \
	--make-bed \
	--memory 8000 \
	--threads 1 \
	--out ${ASSOC}/Var_File_3
	echo \# Pulled Vars for All Genes to 1 file -  `date` \#
	echo `date` "All Genes Pulled to 1 File" >> ${UPDATE_FILE}
fi

###########################################################
## Pull out variants for each individual gene
 # Should be much quicker after making a smaller variant file above
echo \### 4b - `date` \###
echo \### Make BED/012 for each Gene \###
echo `date` "4b - Make BED/012 for each Gene" >> ${UPDATE_FILE}

## Loop thru to pull out variants for each individual gene
 # Output should be .bed file for each gene
for G in `cat ${NEW_GENE_LIST}`
do
	# Make PLINK directory for each gene
	if [ ! -d "${ASSOC}/${G}/PLINK_FILES/" ] ; then mkdir ${ASSOC}/${G}/PLINK_FILES ; fi
	# Output variant PED file for each gene
	${PLINK} \
	--bfile ${ASSOC}/Var_File_3 \
	--silent \
	--hwe 1e-50 midp \
	--extract range ${ASSOC}/${G}/${G}.pl.list \
	--make-bed \
	--out ${ASSOC}/${G}/PLINK_FILES/${G}
	# Output to 012 Format
	${PLINK} \
	--bfile ${ASSOC}/${G}/PLINK_FILES/${G} \
	--silent \
	--recode A \
	--out ${ASSOC}/${G}/PLINK_FILES/${G}
	# Copy to ".012" file location
	cp ${ASSOC}/${G}/PLINK_FILES/${G}.raw ${ASSOC}/${G}/PLINK_FILES/${G}.012
	echo \# 4 - Output .bed for ${G} at `date` \#
	echo `date` "BED/012 Made for: " ${G} >> ${UPDATE_FILE}
done

echo `date` "Done pulling out variants for each Gene" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 5 ## Get Gene Annotations from Cypher/SG-ADVISER File ##
###########################################################
if [ "$START_STEP" -le 5 ]; then
 # Python - Get_All_Gene_Annot.py
echo \### 5 - `date` \###
echo \### Pull of Annotations for each Gene - Tabix \###
echo `date` "5 - Pull out Annotations for each Gene - Tabix" >> ${UPDATE_FILE}

# Use Tabix to Pull out Annotations
zcat ${ANNOTS} | head -1 > ${ASSOC}/Header.txt
for G in `cat ${NEW_GENE_LIST}`
do
	COORDINATES=`cat ${ASSOC}/${G}/${G}.list`
	tabix ${ANNOTS} chr${COORDINATES} > ${ASSOC}/${G}/${G}.Annot.txt
	echo \# 5 - Annotations Pulled for ${G} at `date` \#
	echo `date` "Pulled Annotations for: " ${G} >> ${UPDATE_FILE}
done

echo `date` "Done pulling out Annotations using Tabix" >> ${UPDATE_FILE}
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

Rscript ${GET_FUNC_VARS_R} ${ASSOC} ${NEW_GENE_LIST} # ${ANNOTS}

echo `date` "Done Identifying Functional Variants using Rscript" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 7 ## Make New Covariate File Containing PCs ############
###########################################################
if [ "$START_STEP" -le 7 ]; then
 # If no PCs exist, make them
   # Convert to PED and then calculate PCs
echo \### 7 - `date` \###
echo \### Compile Covariates with PCs \###
echo `date` "7 - Compile Covariates with PCs" >> ${UPDATE_FILE}

# If No Principal Components Exist, Make Them
if [ $EIG_VEC = "F" -a $PC_COUNT -gt 0 -a ${VAR_FILE: -4} == ".vcf" ] ; then
	# Convert Original VCF to PED File
	${VCF_TOOLS} --vcf ${VAR_FILE} --plink --out ${ASSOC}/Full
	# Use PED to run PCA
	${PLINK} \
	--file ${ASSOC}/Full \
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
## 8 ## Run Single-Locus Association on each Variant #####
###########################################################
if [ "$START_STEP" -le 8 ]; then
 # Plink - Logistic/Linear
echo \### 8 - `date` \###
echo \### Running Single-Locus Association - PLINK \###
echo `date` "8 - Running Single-Locus Association - PLINK" >> ${UPDATE_FILE}

for G in `cat ${NEW_GENE_LIST}`
do
	${PLINK} \
	--bfile ${ASSOC}/${G}/PLINK_FILES/${G} \
	--silent \
	--pheno ${PHENO_FILE} \
	--covar ${NEW_COV_FILE} --covar-name ${COVS_COMMAND} \
	--${SUFFIX} hide-covar --adjust \
	--allow-no-sex \
	--maf 0.01 \
	--out ${ASSOC}/${G}/PLINK_FILES/PLINK_${G}
	echo \# Ran Single-Locus Tests for ${G} at `date` \#
	echo `date` "Ran Single-Locus Tests for: " ${G} >> ${UPDATE_FILE}
done

echo `date` "Single-Locus Associations Done" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 9 ## Run SKAT Tests ###################################
###########################################################
if [ "$START_STEP" -le 9 ]; then
 # Rscript
 # Use different collapsed sets
 # Different Tests
echo \### 9 - `date` \###
echo \### Run SKAT Tests - Rscript \###
echo `date` "9 - Run SKAT Tests - Rscript" >> ${UPDATE_FILE}

Rscript ${SKAT_R} ${ASSOC} ${NEW_GENE_LIST} ${PHENO_FILE} ${PHENO_TYPE} ${NEW_COV_FILE} ${COVS_COMMAND} ${GENE_START}

echo `date` "SKAT Tests Done" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 10 ## Make Manhattan Plot for Single Locus Tests ########
###########################################################
if [ "$START_STEP" -le 10 ]; then
 # Rscript
 # Also compiles list of SNPs beyond a given threshold (0.1?)
echo \### 10 - `date` \###
echo \### Make Manhattan Plot of Genes - Rscript \###
echo `date` "10 - Make Manhattan Plot of Genes - Rscript" >> ${UPDATE_FILE}

Rscript ${GENE_MANH} ${NEW_GENE_LIST} ${COV_FILE} ${PHENO_TYPE}

echo `date` "Plotting the Skyline Done" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 11 ## Plot Results of SKAT Tests #######################
###########################################################
if [ "$START_STEP" -le 11 ]; then
 # Rscript
 # TBD...
echo \### 11 - `date` \###
echo \### Plot SKAT Results - Rscript \###
echo `date` "11 - Plot SKAT Results - Rscript" >> ${UPDATE_FILE}

Rscript ${PLOT_SKAT_TESTS_R} ${ASSOC}

echo `date` "Plotting SKAT Results Done" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 12 ## END OF DOC...SO FAR ##############################
###########################################################
echo `date` "DONE!!!" >> ${UPDATE_FILE}
echo \#################################################### >> ${UPDATE_FILE}
echo \#################### DONE!! \######################## >> ${UPDATE_FILE}
printf "DONE!!!"
printf "V\nV\nV\nV\nV\nV\nV\nV\n"









