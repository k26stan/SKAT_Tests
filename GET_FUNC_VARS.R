## Rscript to Compile the Gene Coordinates of multiple genes into a handy table ##
## This table will be used in a python script to get annotations for these genes ##

# Name: GET_FUNC_VARS.R

# Usage:
  # Rscript <This_Script.R> <Path/To/Directory> <Gene_List>
  # Rscript GET_FUNC_VARS.R ${ASSOC} ${GENE_LIST} ${ANNOTS}

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c("/projects/janssen/RV_Rare_Variant/20140401_2013_Plenge_Vars_G_DAS_24_C", "/home/kstandis/MrOSexome/2013_Plenge_Vars.txt")
# LINE <- c("/home/kstandis/MrOSexome/20140403_2012_GWAS_Tb1_Quant_Pheno", "/home/kstandis/MrOSexome/2012_GWAS_Tb1.txt")
# LINE <- c("/projects/ps-jcvi/projects/Kidney/20140506_Test_SKAT_Cap_Reject", "/projects/ps-jcvi/projects/Kidney/20140506_Test_SKAT_Cap_Reject/Alt_Test_SKAT_Cap.txt", "/projects/ps-jcvi/projects/Kidney/SUMMARY_sga-upload_novel_annotation.txt")
PathToBaseDir <- LINE[1]
PathToGeneList <- LINE[2]
PathToAnnots <- LINE[3]

print(PathToBaseDir)
print(PathToGeneList)

########################################################
## Open Gene List and All Gene Coordinates #############
########################################################

## Load Annotation Headers
#HEAD <- as.character(t(read.table(PathToAnnots,sep="\t",nrow=1,header=F))[,1])
HEAD <- as.character(t(read.table(paste(PathToBaseDir,"/Header.txt",sep=""),sep="\t",nrow=1,header=F))[,1])

GENE_LIST <- read.table(PathToGeneList)
print("Loaded Gene List")
#print(dim(GENE_LIST))
NAMES <- as.character(GENE_LIST[,1])

GENES <- list()
for (gene_num in 1:nrow(GENE_LIST)) {
	GENE <- NAMES[gene_num]
	PathToFile <- paste(PathToBaseDir,"/",GENE,"/",GENE,".Annot.txt",sep="")
	if (file.info(PathToFile)$size!=0) {
		GENES[[gene_num]] <- read.table(PathToFile, sep="\t", header=F, fill=T, quote="")
		colnames(GENES[[gene_num]]) <- HEAD
	}
}
names(GENES) <- NAMES
print(summary(GENES))

## Go through Genes, pull out functional Vars, Save Files
print("Looping thru genes!!")
for (gene_num in 1:length(GENES)) {
	DATA <- GENES[[gene_num]]
	GENE <- NAMES[gene_num]
	print(paste(gene_num,GENE,sep=" - "))
	Func_List <- Damg_List <- c()
	gene_col <- as.character(DATA$Gene)
	loc_col <- as.character(DATA$Location)
	cod_col <- as.character(DATA$Coding_Impact)
	for (var in 1:nrow(DATA)) {
		slash_gene <- grep(GENE, strsplit(gene_col[var],"///")[[1]])
		if ( length(grep("Exon",strsplit(loc_col[var],"///")[[1]][slash_gene])) > 0 | length(grep("UTR",strsplit(loc_col[var],"///")[[1]][slash_gene])) > 0 ) {
			#print("Got in the loop!")
			Func_List <- c(Func_List, var)
		}
		if ( length(grep("Nonsynonymous",strsplit(cod_col[var],"///")[[1]][slash_gene])) > 0 | length(grep("Nonsense",strsplit(cod_col[var],"///")[[1]][slash_gene])) > 0 ) {
			#print("Got in the loop!")
			Damg_List <- c(Damg_List, var)
		}
	}
	Func_Coords <- data.frame(gsub("chr","",DATA$Chromosome[Func_List]), DATA$End[Func_List])
	write.table(Func_Coords, paste(PathToBaseDir,"/",GENE,"/",GENE,".Func.list",sep=""), sep="\t", row.names=F, col.names=F, quote=F)
	Damg_Coords <- data.frame(gsub("chr","",DATA$Chromosome[Damg_List]), DATA$End[Damg_List])
	write.table(Damg_Coords, paste(PathToBaseDir,"/",GENE,"/",GENE,".Damg.list",sep=""), sep="\t", row.names=F, col.names=F, quote=F)
}





