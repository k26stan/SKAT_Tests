Here is a description of how to run: Gene_Assoc.sh
  -This file calls several smaller modules to run Rare Variant Analysis on a specified set of Genes

#####################################################
## WHAT YOU'LL NEED TO GET GOING ####################
#####################################################

A) What files/parameters you'll need
  -DATE: This is used in naming the files/directories
    -Helps keep track if you're running multiple times w/ different parameters
  -HOME_DIR: Your test directory will be created in this directory and all subsequent files will be written in the test directory
  -GENE_LIST: A list of Genes that you're interested in
    -If you have coordinates, you can include those as "CHR:START-STOP"
  -IN_VCF: A vcf file with all of your samples and variants
  -ANNOTS: An annotation file from Cypher or SG-Adviser
    -This is used to pull out functional annotations for each gene (e.g., coding regions, non-synonymous variants)
  -PHENO_FILE: A phenotype file in tab-delimited PLINK format
    -(colnames = IID, FID, Pheno)
  -PHENO_TYPE: C for Continuous, B for Binary
  -COV_FILE: (Optional) A file containing all Covariates under consideration
    -Also in PLINK format (colnames = IID, FID, COV_1, COV_2, ..., COV_N)
    -If no Covariates to be included, type F
  -COVS: Covariates you want to include in model. 
    -Will be ignored if no Covariate File is included
  -EIG_VEC: File with Eigenvectors for Variants as calculated/output by Plink -pca
    -If no Eigenvectors are provided, type F
      -PCA will be done and Eigenvectors will be calculated and used as covaraites
  -PC_COUNT: How many Principal Components would you like to include as Covariates?
  -START_STEP: Which Step (see Gene_Assoc.sh) would you like to begin on?
    -This is handy if it fails in the middle for some reason and you need to pick up where you left off
B) What public tools you'll need
  -GATK_JAR: Path to "GenomeAnalysisTK.jar"
  -REF_FA: Path to Reference Genome (same one that your samples were called on)
  -VCF_TOOLS: Path to "vcftools" executable
  -PLINK: Path to "plink" executable
  -GENE_TABLE: A Gene table that I've got housed in my home directory
    -Location: /home/kstandis/HandyStuff/GG-Gene_Names_DB.txt
C) What custom modules will be called by Gene_Assoc.sh
    -(See "Gene_Assoc.sh" for location of files)
    -(Should be in same directory as this README.txt)
  -COMPILE_GENE_COORDS_R: Rscript to Compile Gene Coordinates (KS'14)
  -GET_ALL_GENE_ANNOT_PY: Python script to Pull Annotations from Cypher/SG-ADVISER file(KS'14)
  -GET_FUNC_VARS_R: Rscript to pull identify Functional Variants from Annotations (KS'14)
  -MAKE_COV_TAB_R: Rscript to Make Covariate Table that includes PCs (KS'14)
  -GENE_MANH: Rscript to Make Manhattan Plot from Single-Locus Results for the Specified Genes (KS'14)
  -SKAT_R: Rscript that takes in Variants/Positions and Runs SKAT_O/SKAT_CR tests on a few Collapsed Sets for each Gene (KS'14)
  -PLOT_SKAT_TESTS_R: Rscript that Plots the results of the SKAT tests (KS'14)

#####################################################
## HOW TO EXECUTE Gene_Assoc.sh #####################
#####################################################

1) Set up a QSUB shell with time, jobs names, accounts, etc...
2) In QSUB shell, set up 12 Parameters/Inputs used in the pipeline
  -See 0A (above) for Specifics of Parameters/Inputs
3) Call Gene_Assoc.sh followed by 12 Parameters
  -e.g., /Path/To/Gene_Assoc.sh ${DATE} ${HOME_DIR} ${GENE_LIST} ${IN_VCF} ${ANNOTS} ${PHENO_FILE} ${PHENO_TYPE} ${COV_FILE} ${COVS} ${EIG_VEC} ${PC_COUNT} ${START_STEP}

#####################################################
## HOW Gene_Assoc.sh WORKS ##########################
#####################################################

## Gene_Assoc.sh is broken into separate modules that follow this format:
1) Set Up Paths
  -Take in Specified Parameters
  -Set Paths to Tools, Scripts, and Files
  -Set up other Variables to be used in commands and file/directory names
  -Create a Directory ($ASSOC) where everything will be saved (name includes Gene List, Phenotype, and Covariates)
2) Verify/Retrieve Gene Coordinates
  -Bash
  -Check if Coordinates are provided by GENE_LIST
  -If not, pull them out of GENE_TABLE for each Gene
  -Create directories within the $ASSOC directory for each Gene and make a list for coordinates within each Gene Directory
3) Compile Gene Coordinates
  -Rscript
  -Compiles Gene Coordinates into a Handy Table and gets rid of Genes where coordinates are unknown
  -Also adds 10K bp to each end of the coordinates
4) Pull out SNPs for each Gene
  -GATK
  -Save a separate VCF for each Gene (within it's own directory)
5) Get Gene Annotatiosn from Cypher/SG-ADVISER File
  -Python
  -Pull out Annotations for the specified coordinates for each Gene
6) Pull out Functional Variants from Annotations
  -Rscript
  -Go thru each Gene, load the Annotations, and Pull out Variants that fit various criteria
  -Func.list: UTR or Exonic
  -Damg.list: Non-synonymous Protein Coding Variants
7) Re-Genotype Missing Variants to 0/0
  -Bash
  -As a result of Merging several VCF files, some positions are considered missing from the VCF file
  -I use "sed" to replace "./." with "0/0" as we assume that it's absence from the VCF file suggests it's homozygous reference
  -Leaving it as "missing" would cause problems later on when population allele frequencies are considered (PCA, Association)
8) Convert Files to Plink Format
  -VCFtools
  -Convert VCF file to PED/MAP files
  -Also converts to 0/1/2 format (used in SKAT tests)
9) Make New Covariate File Containing PCs
  -Plink (If no PCs are provided, they will be calculated using Plink --pca)
  -Rscript
  -Compiles previously specified Covariate file with PCs
  -If no Covariates were previously indicated, PC file will be used as covariate file
10)Run Single-Locus Association on each Variant
  -Plink
  -Phenotype considered was one of the original files specified
  -Includes Clinical Covariates and PCs in Linear/Logistic Regression for each individual SNP
  -MAF is specified at 1% (Anything w/ MAF<1% will not be tested)
11)Make Manhattan Plot for Single-Locus Tests
  -Rscript
  -Takes in output files from Plink and spits out a single Manhattan and QQ plot with P-values for all the Genes
  -Also outputs a list for each gene with variants meeting a specific p-value threshold (p<0.1)
    -To be used in collapsed SKAT tests
12)Run SKAT Tests
  -Rscript
  -Takes in .012 files for each gene, creates a few collapsed sets, performs SKAT test
    -Full: All Rare Variants within Gene Limits + 10Kb on either side
    -Edge: All Rare Variants within Gene Limits
    -Func: Rare Variants in Exons and UTRs
    -Damg: Non-Synonymous Rare Variants
    -CRFull: All Common and Rare Variants within Gene Limits + 10Kb on either side
    -Thrsh: Common and Rare Variants that 1) have p<0.1 in Single-Locus test [or] 2) Are in "Damg" list
  -Outputs a single, compiled SKATout.Rdata object to be loaded/used when plotting these (next step)
13)Plots Results of SKAT Tests
  -Rscript
  -Outputs plots for the p-values, number of variants included in each test, and a QQ plot for gene-level tests

#####################################################
## EXAMPLE SUBMISSION ###############################
#####################################################

## Names/Paths
echo \### Defining Set Variables and Paths at `date` \###
DATE=20140614
HOME_DIR=/projects/janssen/RV_Rare_Variant

## Files
GENE_LIST=${HOME_DIR}/2014_Plenge.txt
IN_VCF=/projects/janssen/VCFs/5-MERGE/HC_FULL.SNP.vcf
ANNOTS=/projects/janssen/ANNOTATE/JnJ_121613_all_annotations.txt.gz
PHENO_FILE=/projects/janssen/ASSOCIATION/PH-PHENOTYPES/LT8_DEL_28_BL.txt
PHENO_TYPE=C
COV_FILE=/projects/janssen/ASSOCIATION/PH-PHENOTYPES/COV.txt
COVS=`echo DAS_BL AGE SEX`
EIG_VEC=/projects/janssen/ASSOCIATION/EIGEN/HC_FULL.eigenvec
PC_COUNT=2
START_STEP=1

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

