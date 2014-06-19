# Make Covariate File w/ EigenVectors #

# Make_Cov_Tab.R <Path/To/EigenVectors> <Path/To/Covariate_Table> <Set>

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c("/projects/janssen/scripts/AS-ASSOCIATION/20140514pca_HT_SEX_AGE/20140514pca-BGI_HT_.eigenvec","/projects/janssen/scripts/AS-ASSOCIATION/COV_Tables/All_BigCall.txt","BGI")
PathToVec <- LINE[1]
PathToCov <- LINE[2]
PathToAssoc <- LINE[3]


## Load Tables
VEC <- read.table(PathToVec,header=T)
COV <- read.table(PathToCov,sep="\t",header=T)

## Merge Tables
VEC <- VEC[,c(1,3:ncol(VEC))]
MRG <- merge(x=COV,y=VEC,by="FID")

## Write Table
write.table(MRG,PathToAssoc,sep="\t",row.names=F,col.names=T,quote=F)

######### END OF DOC ###########
