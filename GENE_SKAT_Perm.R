## Rscript to run CEP-SKAT on Variants ##
## This one will be specific to MrOSexome ##

# Name: GENE_SKAT.R

# Usage:
  # Rscript ${SKAT_R} ${ASSOC} ${NEW_GENE_LIST} ${PHENO_FILE} ${PHENO_TYPE} ${NEW_COV_FILE} ${COVS_COMMAND}

#############################################
## PARSE COMMAND LINE #######################
#############################################

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c("/Users/kstandis/Data/MrOSexome/20140612_TEST_GENES_Quant_Pheno","/Users/kstandis/Data/MrOSexome/20140612_TEST_GENES_Quant_Pheno/Alt_TEST_GENES.txt","/Users/kstandis/Data/MrOSexome/DATA/Quant_Pheno.txt","C","Data/MrOSexome/20140612_TEST_GENES_Quant_Pheno/Cov_w_PCs.txt","PC1,PC2,PC3")
# LINE <- c("/home/kstandis/MrOSexome/20140613b_2012_GWAS_Tb2_Quant_Pheno","/home/kstandis/MrOSexome/20140613b_2012_GWAS_Tb2_Quant_Pheno/Alt_2012_GWAS_Tb2.txt", "/home/kstandis/MrOSexome/PHENOTYPE/Quant_Pheno.txt","C", "/home/kstandis/MrOSexome/20140613b_2012_GWAS_Tb2_Quant_Pheno/Cov_w_PCs.txt","PC1,PC2,PC3")
# LINE <- c("/projects/janssen/RV_Rare_Variant/20140616_2013_Plenge_Vars_LT8_DEL_28_BL","/projects/janssen/RV_Rare_Variant/20140616_2013_Plenge_Vars_LT8_DEL_28_BL/Alt_2013_Plenge_Vars.txt","/projects/janssen/RV_Rare_Variant/../ASSOCIATION/PH-PHENOTYPES/LT8_DEL_28_BL.txt","C","/projects/janssen/RV_Rare_Variant/20140616_2013_Plenge_Vars_LT8_DEL_28_BL/Cov_w_PCs.txt","DAS_BL,AGE,SEX,PC1,PC2")
PATH <- LINE[1] # Path To General Directory
G_LIST <- LINE[2] # Path to Gene List
PHE <- LINE[3] # Path to Phenotype File
PHENO_TYPE <- LINE[4] # (B)inary or (C)ontinuous Phenotype?
COV <- LINE[5] # Path to Covariate File
COVS <- LINE[6] # List of Covariates to Include
print("############# Command Line Parsed ###############")
print(paste("Path to Base Directory is:",PATH))
print(paste("Gene List is:",G_LIST))
print(paste("Path To Pheno File is:",PHE))
print(paste("Phenotype is:",PHENO_TYPE))
print(paste("Path to Covariates is:",COV))
print(paste("Including Covariates:",COVS))

## Split Arguments to get out relevant info
Dir <- strsplit(PATH,"/")
Dir <- Dir[[1]][grep("2014",Dir[[1]])]
Date <- strsplit(Dir,"_")[[1]][1]
Home <- paste(strsplit(PATH,"/")[[1]][2:4],collapse="/")
PathToUpdate <- paste(PATH,"/Update.txt",sep="")

## Set PHENO_TYPE Variable (binary or continuous) 
if (PHENO_TYPE == "B") { PHENO_TYPE <- "D" }

#############################################
## LOAD PHENO/COV DATA ######################
#############################################
## Load SKAT Libraries
library(SKAT)
library(CEPSKAT)

## Load Phenotype Files
print("Loading Phenotype Information")
P <- read.table(PHE, header=T, colClasses=c("character","character","numeric"))
print("Phenotype File Loaded")

## Load Covariate Files
X_Arr <- read.table(COV,header=T)
print("Covariate File Loaded")

## Get relevant Sample Names
Samp_Names <- intersect(P$IID[which(P$Pheno!=-9)],X_Arr$IID)

## Filter and Re-format Pheno/Covar Files
P <- P[which(P$IID %in% Samp_Names),] ; P <- P[order(P$IID),]
Cov_Names <- strsplit(COVS,",")[[1]]
X_Arr <- X_Arr[which(X_Arr$IID %in% Samp_Names),] ; X_Arr <- X_Arr[order(X_Arr$IID),]
X <- data.matrix(X_Arr[,Cov_Names]) ; rownames(X) <- X_Arr$IID
print("Pheno/Cov Data Reformatted")

#############################################
## CALCULATE NULL MODELS ####################
#############################################

## Null model for covariates
Null <- SKAT_Null_Model(P$Pheno ~ X, out_type=PHENO_TYPE) # Choose whichever covariates I want
print("Null Models Calculated")

##############################################################
## BEGIN LOOPING THRU GENES ##################################
##############################################################
## For each Gene, I want:
 # SKAT Test Results
   # SKAT_O - Full,Full+10K,Func,Damg ("Full" has +10K, "Edge" is at the Edge of the Gene)
   # SKAT_RC - Thrsh (includes 10K)
 # Permuted P-values for each Test Above
   # 100x to start? I dunno...just set it up then talk to Nik
 # Number of Variants
   # Full, Full+10K, Func, Damg, Thrsh
 # Rho Values for each Actual Test
 # TBD...

## Load Gene List and Coordinates
GENE_TAB <- as.character(read.table(G_LIST)[,1])
GENE_COORD <- read.table(paste(PATH,"/Gene_Coords.txt",sep=""),sep="\t",header=T)
print("Gene List and Coordinates Loaded")

## Number of Permutations
PERMS <- 100 # 100 takes ~18+ seconds for 100 Vars...give or take...maybe?

## Set up Variable to Fill in whilst looping thru the Genes
 # For SKAT-O P-Values
SO_Full <- SO_Edge <- SO_Func <- SO_Damg <- SO_CRFull <- SO_Thrsh <- rep(NA,length(GENE_TAB))
names(SO_Full) <- names(SO_Edge) <- names(SO_Func) <- names(SO_Damg) <- names(SO_CRFull) <- names(SO_Thrsh) <- GENE_TAB
 # For Strictly SKAT P-Values (Rho=0)
SK_Full <- SK_Edge <- SK_Func <- SK_Damg <- rep(NA,length(GENE_TAB))
names(SK_Full) <- names(SK_Edge) <- names(SK_Func) <- names(SK_Damg) <- GENE_TAB
 # For Burden Test P-Values (Rho=1)
BD_Full <- BD_Edge <- BD_Func <- BD_Damg <- rep(NA,length(GENE_TAB))
names(BD_Full) <- names(BD_Edge) <- names(BD_Func) <- names(BD_Damg) <- GENE_TAB
 # For Permuted P-Values
PM_Full <- PM_Edge <- PM_Func <- PM_Damg <- PM_CRFull <- PM_Thrsh <- array(,dim=c(length(GENE_TAB),PERMS))
rownames(PM_Full) <- rownames(PM_Edge) <- rownames(PM_Func) <- rownames(PM_Damg) <- rownames(PM_CRFull) <- rownames(PM_Thrsh) <- GENE_TAB
colnames(PM_Full) <- colnames(PM_Edge) <- colnames(PM_Func) <- colnames(PM_Damg) <- colnames(PM_CRFull) <- colnames(PM_Thrsh) <- paste("Perm",1:PERMS,sep="")
 # For Number of Variants
NV_Full <- NV_Edge <- NV_Func <- NV_Damg <- rep(NA,length(GENE_TAB))
names(NV_Full) <- names(NV_Edge) <- names(NV_Func) <- names(NV_Damg) <- GENE_TAB
NV_CRFull <- NV_Thrsh <- array(,dim=c(length(GENE_TAB),2))
rownames(NV_CRFull) <- rownames(NV_Thrsh) <- GENE_TAB
colnames(NV_CRFull) <- colnames(NV_Thrsh) <- c("Common","Rare")
 # For Rho Values
RO_Full <- RO_Edge <- RO_Func <- RO_Damg <- RO_CRFull <- RO_Thrsh <- rep(NA,length(GENE_TAB))
names(RO_Full) <- names(RO_Edge) <- names(RO_Func) <- names(RO_Damg) <- names(RO_CRFull) <- names(RO_Thrsh) <- GENE_TAB

start_gene_loop <- proc.time()
## Loop Through Genes and do this Shiz...
for (i in 1:length(GENE_TAB)) {
	print(paste("Running Skat on Gene",i,"of",length(GENE_TAB)))
	GENE <- GENE_TAB[i]
	write(paste(Sys.time(),"- Running Skat on:",GENE,"-",i,"of",length(GENE_TAB)),PathToUpdate,append=T)

#############################################
## SET PATHS and LOAD GENOTYPE DATA #########
#############################################

## Set Paths to Genotype and Position Files
Full_Array <- paste(PATH, "/", GENE, "/", GENE, ".012", sep="")
Func_Pos_Nm <- paste(PATH, "/", GENE, "/", GENE, ".Func.list", sep="")
Damg_Pos_Nm <- paste(PATH, "/", GENE, "/", GENE, ".Damg.list", sep="")
Thrsh_Pos_Nm <- paste(PATH, "/", GENE, "/", GENE, ".Thrsh.list", sep="")

## Load Sample and Variant ID Files
print(paste("Gene:",GENE,"-","Loading Genotype Information"))
G_Names <- as.character(read.table(paste(Full_Array,".indv",sep=""), sep="\t", header=F)[,1])
G_Locs <- read.table(paste(Full_Array,".pos",sep=""), sep="\t", header=F, colClasses=c("factor","numeric"))

File_Flag <- 0
if (file.exists(Full_Array)==T) {
	File_Flag <- 1
	G_Full <- read.table(Full_Array, sep="\t", header=F)
	G_Full <- G_Full[,2:ncol(G_Full)]
	colnames(G_Full) <- G_Locs[,2] #paste(G_Locs[,1],G_Locs[,2],sep=":")
	rownames(G_Full) <- G_Names
	G_Full <- G_Full[order(rownames(G_Full)),] }
if (file.exists(Func_Pos_Nm)==T) {
	if (length(readLines(Damg_Pos_Nm))>0) {
		Func_Locs <- read.table(Func_Pos_Nm, sep="\t",header=F)
		File_Flag <- c(File_Flag,2)
	}
}
if (file.exists(Damg_Pos_Nm)==T) {
	if (length(readLines(Damg_Pos_Nm))>0) {
		Damg_Locs <- read.table(Damg_Pos_Nm, sep="\t",header=F)
		File_Flag <- c(File_Flag,3)
	}
}
if (file.exists(Thrsh_Pos_Nm)==T) {
	if (length(readLines(Damg_Pos_Nm))>0) {
		Thrsh_Locs <- read.table(Thrsh_Pos_Nm, sep="\t",header=F)
		File_Flag <- c(File_Flag,4)
	}
}

## Check to Make Sure all the Sample names still Match
if (length(which(sort(G_Names)!=sort(Samp_Names)))!=0 ) { Samp_Names <- intersect(G_Names,Samp_Names) }

## Make 1 the Minor-Allele
G_Full_Flip <- which( apply(G_Full,2,sum)/(2*nrow(G_Full)) > .5 )
if (length(G_Full_Flip)>0) {
	G_Full[,G_Full_Flip] <- -G_Full[,G_Full_Flip]+2
	rownames(G_Full) <- G_Names
}

print("BEGINNING SKAT TESTS")
#############################################
## SKAT-O for FULL DATA #####################
#############################################
if (1 %in% File_Flag) {
	G1 <- as.matrix(G_Full[Samp_Names,])
	# Pull out only Rare Variants
	AF <- apply(G1,2,sum)/(2*nrow(G1))
	G_Rare <- G1[,which( AF <= 1/sqrt(2*nrow(G1)) & AF > 0 )]
	print(paste("Gene:",GENE,"-",ncol(data.frame(G_Rare)),"Rare Variants"))
	# Run SKAT-O Test
	if (ncol(data.frame(G_Rare))>1) {
		MOD <- SKAT(G_Rare, Null, kernel="linear.weighted", method="optimal.adj", weights.beta=c(1,25))
		# Compile Results (SO/SK/BD/PM/NV/RO)
		SO_Full[i] <- MOD$p.value
		SK_Full[i] <- MOD$param$p.val.each[which(MOD$param$rho==0)]
		BD_Full[i] <- MOD$param$p.val.each[which(MOD$param$rho==1)]
		NV_Full[i] <- MOD$param$n.marker.test
		RO_Full[i] <- MOD$param$rho_est
		start <- proc.time()
		print(paste("Gene:",GENE,"-","Permuting Full Set"))
		for (j in 1:PERMS) {
			G_Perm <- G_Rare[sample(1:nrow(G_Rare),nrow(G_Rare),replace=F),]
			P_MOD <- SKAT(G_Perm, Null, kernel="linear.weighted", method="optimal.adj", weights.beta=c(1,25))
			PM_Full[i,j] <- P_MOD$p.value
		}
		proc.time()-start
	}

#############################################
## SKAT-O for EDGE DATA #####################
#############################################
	MIN <- GENE_COORD$Min[which(GENE_COORD$Gene==GENE)]+10000
	MAX <- GENE_COORD$Max[which(GENE_COORD$Gene==GENE)]-10000
	WHICH <- which(as.numeric(colnames(G_Rare))>=MIN & as.numeric(colnames(G_Rare))<=MAX)
	if (!identical(WHICH,1:ncol(G_Rare))) { # If there are variants in the 10K region on either end of the Gene
		G <- as.matrix(G_Rare[Samp_Names,WHICH])
		print(paste("Gene:",GENE,"-",ncol(data.frame(G)),"Rare Variants in Edge Set"))
		# Run SKAT-O Test
		if (ncol(data.frame(G))>1) {
			MOD <- SKAT(G, Null, kernel="linear.weighted", method="optimal.adj", weights.beta=c(1,25))
			# Compile Results (SO/SK/BD/PM/NV/RO)
			SO_Edge[i] <- MOD$p.value
			SK_Edge[i] <- MOD$param$p.val.each[which(MOD$param$rho==0)]
			BD_Edge[i] <- MOD$param$p.val.each[which(MOD$param$rho==1)]
			NV_Edge[i] <- MOD$param$n.marker.test
			RO_Edge[i] <- MOD$param$rho_est
			start <- proc.time()
			print(paste("Gene:",GENE,"-","Permuting Edge Set"))
			for (j in 1:PERMS) {
				G_Perm <- G[sample(1:nrow(G),nrow(G),replace=F),]
				P_MOD <- SKAT(G_Perm, Null, kernel="linear.weighted", method="optimal.adj", weights.beta=c(1,25))
				PM_Full[i,j] <- P_MOD$p.value
			}
			proc.time()-start
		}
	}
	else { # If there were no variants within 10K of the edge.
		# Take Values from the "Full" version
		print(paste("Gene:",GENE,"-","Edge Set the same as Full Set"))
		SO_Edge[i] <- SO_Full[i]
		SK_Edge[i] <- SK_Full[i]
		BD_Edge[i] <- BD_Full[i]
		NV_Edge[i] <- NV_Full[i]
		RO_Edge[i] <- RO_Full[i]
		PM_Edge[i,] <- PM_Full[i,]
	}
} # Closes: if (1 %in% File_Flag)

#############################################
## SKAT-O for FUNC DATA #####################
#############################################
if (2 %in% File_Flag) {
	WHICH <- which(colnames(G_Rare) %in% Func_Locs[,2])
	G <- as.matrix(G_Rare[Samp_Names,WHICH])
	print(paste("Gene:",GENE,"-",ncol(data.frame(G)),"Rare Variants in Func Set"))
	# Run SKAT-O Test
	if (ncol(data.frame(G))>1) {
		MOD <- SKAT(G, Null, kernel="linear.weighted", method="optimal.adj", weights.beta=c(1,25))
		# Compile Results (SO/SK/BD/PM/NV/RO)
		SO_Func[i] <- MOD$p.value
		SK_Func[i] <- MOD$param$p.val.each[which(MOD$param$rho==0)]
		BD_Func[i] <- MOD$param$p.val.each[which(MOD$param$rho==1)]
		NV_Func[i] <- MOD$param$n.marker.test
		RO_Func[i] <- MOD$param$rho_est
		start <- proc.time()
		print(paste("Gene:",GENE,"-","Permuting Func Set"))
		for (j in 1:PERMS) {
			G_Perm <- G[sample(1:nrow(G),nrow(G),replace=F),]
			P_MOD <- SKAT(G_Perm, Null, kernel="linear.weighted", method="optimal.adj", weights.beta=c(1,25))
			PM_Func[i,j] <- P_MOD$p.value
		}
		proc.time()-start
	}
} # Closes: if (2 %in% File_Flag)

#############################################
## SKAT-O for DAMG DATA #####################
#############################################
if (3 %in% File_Flag) {
	WHICH <- which(colnames(G_Rare) %in% Damg_Locs[,2])
	G <- as.matrix(G_Rare[Samp_Names,WHICH])
	print(paste("Gene:",GENE,"-",ncol(data.frame(G)),"Rare Variants in Damg Set"))
	# Run SKAT-O Test
	if (ncol(data.frame(G))>1) {
		MOD <- SKAT(G, Null, kernel="linear.weighted", method="optimal.adj", weights.beta=c(1,25))
		# Compile Results (SO/SK/BD/PM/NV/RO)
		SO_Damg[i] <- MOD$p.value
		SK_Damg[i] <- MOD$param$p.val.each[which(MOD$param$rho==0)]
		BD_Damg[i] <- MOD$param$p.val.each[which(MOD$param$rho==1)]
		NV_Damg[i] <- MOD$param$n.marker.test
		RO_Damg[i] <- MOD$param$rho_est
		start <- proc.time()
		print(paste("Gene:",GENE,"-","Permuting Damg Set"))
		for (j in 1:PERMS) {
			G_Perm <- G[sample(1:nrow(G),nrow(G),replace=F),]
			P_MOD <- SKAT(G_Perm, Null, kernel="linear.weighted", method="optimal.adj", weights.beta=c(1,25))
			PM_Damg[i,j] <- P_MOD$p.value
		}
		proc.time()-start
	}
} # Closes: if (3 %in% File_Flag)

#############################################
## SKAT-CR for FULL DATA ####################
#############################################
if (1 %in% File_Flag) {
	G <- as.matrix(G_Full[Samp_Names,])
	AF <- apply(G,2,sum)/(2*nrow(G))
	G_Com <- G[,which( AF > 1/sqrt(2*nrow(G)) & AF > 0 )]
	print(paste("Gene:",GENE,"-",ncol(data.frame(G_Com)),"Common Variants in Full Set"))
	if (ncol(data.frame(G))>1) {
		# Get Rho Values for Common/Rare Variants
		rho_R <- RO_Full[i]
		rho_C <- 0
		if (ncol(data.frame(G_Com))>1) {
			MOD_O_Com <- SKAT(G_Com, Null, kernel="linear.weighted", method="optimal.adj", weights.beta=c(.5,.5))
			rho_C <- MOD_O_Com$param$rho_est
		}
		# Run SKAT-Common/Rare
		MOD <- SKAT_CommonRare(G, Null, weights.beta.rare=c(1,25), weights.beta.common=c(.5,.5), method="C", r.corr.rare=rho_R, r.corr.common=rho_C)
		# Compile Results (SO/SK/BD/PM/NV/RO)
		SO_CRFull[i] <- MOD$p.value
		NV_CRFull[i,1] <- MOD$n.common
		NV_CRFull[i,2] <- MOD$n.rare
		RO_CRFull[i] <- rho_C
		start <- proc.time()
		print(paste("Gene:",GENE,"-","Permuting C/R Full Set"))
		for (j in 1:PERMS) {
			G_Perm <- G[sample(1:nrow(G),nrow(G),replace=F),]
			P_MOD <- SKAT_CommonRare(G_Perm, Null, weights.beta.rare=c(1,25), weights.beta.common=c(.5,.5), method="C", r.corr.rare=rho_R, r.corr.common=rho_C)
			PM_CRFull[i,j] <- P_MOD$p.value
		}
		proc.time()-start
	}
} # Closes: if (1 %in% File_Flag)

#############################################
## SKAT-CR for THRSH DATA ###################
#############################################
if (4 %in% File_Flag) {
	WHICH <- which(colnames(G_Full) %in% Thrsh_Locs[,2])
	G <- as.matrix(G_Full[Samp_Names,WHICH])
	AF <- apply(G,2,sum)/(2*nrow(G))
	G_Rare <- G[,which( AF <= 1/sqrt(2*nrow(G)) & AF > 0 )]
	G_Com <- G[,which( AF > 1/sqrt(2*nrow(G)) & AF > 0 )]
	print(paste("Gene:",GENE,"-",ncol(data.frame(G_Rare)),"Rare Variants in Thrsh Set"))
	print(paste("Gene:",GENE,"-",ncol(data.frame(G_Com)),"Common Variants in Thrsh Set"))
	if (ncol(data.frame(G))>1) {
		# Get Rho Values for Common/Rare Variants
		rho_R <- 0
		if (ncol(data.frame(G_Rare))>1) {
			MOD_O_Rare <- SKAT(G_Rare, Null, kernel="linear.weighted", method="optimal.adj", weights.beta=c(1,25))
			rho_R <- MOD_O_Rare$param$rho_est
		}
		rho_C <- 0
		if (ncol(data.frame(G_Com))>1) {
			MOD_O_Com <- SKAT(G_Com, Null, kernel="linear.weighted", method="optimal.adj", weights.beta=c(.5,.5))
			rho_C <- MOD_O_Com$param$rho_est
		}
		# Run SKAT-Common/Rare
		MOD <- SKAT_CommonRare(G, Null, weights.beta.rare=c(1,25), weights.beta.common=c(.5,.5), method="C", r.corr.rare=rho_R, r.corr.common=rho_C)
		# Compile Results (SO/SK/BD/PM/NV/RO)
		SO_Thrsh[i] <- MOD$p.value
		NV_Thrsh[i,1] <- MOD$n.common
		NV_Thrsh[i,2] <- MOD$n.rare
		RO_Thrsh[i] <- rho_C
		start <- proc.time()
		print(paste("Gene:",GENE,"-","Permuting C/R Thrsh Set"))
		for (j in 1:PERMS) {
			G_Perm <- G[sample(1:nrow(G),nrow(G),replace=F),]
			P_MOD <- SKAT_CommonRare(G_Perm, Null, weights.beta.rare=c(1,25), weights.beta.common=c(.5,.5), method="C", r.corr.rare=rho_R, r.corr.common=rho_C)
			PM_Thrsh[i,j] <- P_MOD$p.value
		}
		proc.time()-start
	}
} # Closes: if (4 %in% File_Flag)

##############################################################
## CLOSE GENE LOOP ###########################################
##############################################################
print(paste("Gene:",GENE,"-","Elapsed Time:"))
print(proc.time()-start_gene_loop)
} # Closes Gene Loop
print("SKAT TESTS DONE!!")

#############################################
## COMPILE AND SAVE DATA ####################
#############################################
print("COMPILING AND SAVING RESULTS!!")

## Give Names to each Collapsed Set
SET_NAMES <- c("Full","Edge","Func","Damg","CRFull","Thrsh")

## Compile Data into Single Variable (SO/SK/BD/PM/NV/RO)
P_SO <- list(SO_Full,SO_Edge,SO_Func,SO_Damg,SO_CRFull,SO_Thrsh)
P_SK <- list(SK_Full,SK_Edge,SK_Func,SK_Damg)
P_BD <- list(BD_Full,BD_Edge,BD_Func,BD_Damg)
P_PM <- list(PM_Full,PM_Edge,PM_Func,PM_Damg,PM_CRFull,PM_Thrsh)
N_VAR <- list(NV_Full,NV_Edge,NV_Func,NV_Damg,NV_CRFull,NV_Thrsh)
RO_VAL <- list(RO_Full,RO_Edge,RO_Func,RO_Damg,RO_CRFull,RO_Thrsh)
names(P_SO) <- names(P_PM) <- names(N_VAR) <- names(RO_VAL) <- SET_NAMES
names(P_SK) <- names(P_BD) <- SET_NAMES[1:4]
COMPILE <- list(P_SO,P_SK,P_BD,P_PM,N_VAR,RO_VAL)
names(COMPILE) <- c("P_SKAT-O","P_SKAT","P_Burden","P_Permuted","Num_Vars","Rho_Val")

## Save the Compiled data as an object in the overall directory
print(paste("Saving File to Disc: ",PATH,"/SKATout.Rdata",sep=""))
save(COMPILE,file=paste(PATH, "/SKATout.Rdata", sep=""))

###############################################################
#######################    END OF DOC    ######################
###############################################################