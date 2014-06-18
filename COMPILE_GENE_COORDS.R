## Rscript to Compile the Gene Coordinates of multiple genes into a handy table ##
## This table will be used in a python script to get annotations for these genes ##

# Name: COMPILE_GENE_COORDS.R

# Usage:
  # Rscript <This_Script.R> <Path/To/Directory> <Gene_List>
  # Rscript COMPILE_GENE_COORDS.R ${ASSOC} ${GENE_LIST}

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c("/projects/janssen/RV_Rare_Variant/20140401_2013_Plenge_Vars_G_DAS_24_C", "/home/kstandis/MrOSexome/2013_Plenge_Vars.txt")
# LINE <- c("/home/kstandis/MrOSexome/20140403_2012_GWAS_Tb1_Quant_Pheno", "/home/kstandis/MrOSexome/2012_GWAS_Tb1.txt")
# LINE <- c("/home/kstandis/MrOSexome/20140612_TEST_GENES_Quant_Pheno", "/home/kstandis/MrOSexome/20140612_TEST_GENES_Quant_Pheno/TEST_GENES.txt")
PathToBaseDir <- LINE[1]
PathToGeneList <- LINE[2]

print(PathToBaseDir)
print(PathToGeneList)

########################################################
## Open Gene List and All Gene Coordinates #############
########################################################

GENE_LIST <- read.table(PathToGeneList)
print("Loaded Gene List")
print(dim(GENE_LIST))

GENES <- list()
for (gene_num in 1:nrow(GENE_LIST)) {
	GENE <- as.character(GENE_LIST[gene_num,1])
	PathToFile <- paste(PathToBaseDir,"/",GENE,"/",GENE,".list",sep="")
	if (file.info(PathToFile)$size!=0) {
		GENES[[gene_num]] <- read.table(PathToFile)
	}
}
#print(summary(GENES))

########################################################
## Compile Coordinates for each Gene ###################
########################################################

COMPILE_1 <- array( ,dim=c(length(GENES),4) )

for (gene_num in 1:length(GENES)) {
	if (length(dim(GENES[[gene_num]]))) {
		GENE <- as.character(GENE_LIST[gene_num,1])
		DATA <- as.character(GENES[[gene_num]][,1])
		Split1 <- t(sapply(strsplit(DATA,":"),"[",c(1:2)))
		CHROMs <- Split1[,1]
		Split2 <- t(sapply(strsplit(Split1[,2],"-"),"[",c(1:2)))
		MIN <- min(as.numeric(Split2[,1])) - 10000
		MAX <- max(as.numeric(Split2[,2])) + 10000
		if ( length(which(duplicated(CHROMs)))==(length(CHROMs)-1) ) {
			CHROM <- CHROMs[1]
			COMPILE_1[gene_num,] <- c(GENE,CHROM,MIN,MAX)
		}
	}
}
#print(COMPILE_1)

# Re-order
COMPILE_noX <- COMPILE_1[which(COMPILE_1[,2] %in% 1:22),]
COMPILE_XY <- COMPILE_1[which(COMPILE_1[,2] %in% c("X","Y","M")),]

COMPILE_noX <- COMPILE_noX[order(as.character(COMPILE_noX[,2]),as.numeric(COMPILE_noX[,3])),]
COMPILE <- rbind(COMPILE_noX, COMPILE_XY)
colnames(COMPILE) <- c("Gene", "Chrom", "Min", "Max")
print(COMPILE)

# Save this tables to be used in Python script
PathToSave_1 <- paste(PathToBaseDir,"/Gene_Coords.txt",sep="")
print(PathToSave_1)
write.table(COMPILE,PathToSave_1,sep="\t",row.names=F,col.names=T,quote=F)
print("Hopefully that saved...")

# Save new Gene List w/o empty Genes
Split <- strsplit(PathToGeneList,"/")
Gene_List <- Split[[1]][length(Split[[1]])]
NEW_GENE_LIST <- COMPILE[,1]
PathToSave_2 <- paste(PathToBaseDir,"/Alt_",Gene_List,sep="")
print(PathToSave_2)
write.table(NEW_GENE_LIST,PathToSave_2,row.names=F,col.names=F,quote=F)
print("Hopefully the last one saved too...")

# Re-write "GENE.list" files to be used in GATK
for (gene_num in 1:length(NEW_GENE_LIST)) {
	GENE <- NEW_GENE_LIST[gene_num]
	PathToFile <- paste(PathToBaseDir,"/",GENE,"/",GENE,".list",sep="")
	COORDS <- paste(COMPILE[gene_num,3:4],collapse="-")
	OUTPUT <- paste(COMPILE[gene_num,2],COORDS,sep=":")
	write.table(OUTPUT,PathToFile,sep="\t",row.names=F,col.names=F,quote=F)
}

