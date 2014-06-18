## Rscript to Parse through the Rdata output files of SKAT test and make plots ##

# Name: PLOT_SKAT_TESTS.R

# Usage:
  # Rscript <This_Script.R> <Path/To/Directory> <Gene_List>
  # Rscript PLOT_SKAT_TESTS.R ${ASSOC} ${GENE_LIST}

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c("/projects/janssen/RV_Rare_Variant/20140429_2013_Plenge_Vars_B_CAT_52", "/projects/janssen/RV_Rare_Variant/20140429_2013_Plenge_Vars_B_CAT_52/Alt_2013_Plenge_Vars.txt")
# LINE <- c("/home/kstandis/MrOSexome/20140421_2012_GWAS_Unq_Quant_Pheno", "/home/kstandis/MrOSexome/2012_GWAS_Unq.txt")
# LINE <- c("/home/kstandis/MrOSexome/20140613b_2012_GWAS_Tb2_Quant_Pheno")
PathToBaseDir <- LINE[1]

########################################################
## Load SKAT Results ###################################
########################################################

## Path to Skat Results
PathToResults <- paste(PathToBaseDir,"/SKATout.Rdata",sep="")

load(PathToResults)

## Summary of Results that were just Loaded ##
## COMPILE = List of 6 items
 # [[1]]P_SKAT-O = List of P-Values for SKAT-O test for 6 Collapsed Sets
   # [[1]]Full = All Rare Variants +10K on either end of Gene Border
   # [[2]]Edge = All Rare Variants within Gene Borders (w/o extra 10K)
   # [[3]]Func = Rare Exonic and UTR Variants
   # [[4]]Damg = Rare Nonsense and Missense Exonic Variants
   # [[5]]CRFull = Common/Rare SKAT on Full Set
   # [[6]]Thrsh = Common/Rare SKAT on union of "Damg" variants and any Variants with individual association of p<=0.1
 # [[2]]P_SKAT = List of P-Values for SKAT test for 4 Collapsed Sets
   # [[1]]Full = All Rare Variants +10K on either end of Gene Border
   # [[2]]Edge = All Rare Variants within Gene Borders (w/o extra 10K)
   # [[3]]Func = Rare Exonic and UTR Variants
   # [[4]]Damg = Rare Nonsense and Missense Exonic Variants
 # [[3]]P_Burden = List of P-Values for Burden test
   # 4 Collapsed Sets
 # [[4]]P_SKAT-O = Array of Permuted P-Values for SKAT-O test
   # 6 Collapsed Sets
 # [[5]]Num_Vars = Number of Variants included in SKAT-O test
   # Vector for 4 Collapsed Sets
   # Array for 2 Collapsed Sets (counts of Common and Rare Variants included)
 # [[6]]Rho_Val = Optimal Rho Value for SKAT-O test
   # Values for Rare Variant SKAT-O Test for 4 collapsed Sets
   # Values for Common Variant SKAT-O Test for 2 collapsed Sets

## Designate Colors for each Collapsed Set
SET_COLS <- c("chartreuse1","dodgerblue1","goldenrod1","firebrick1","purple1","chocolate1")
# barplot(t(array(1:12,dim=c(6,2))),col=SET_COLS, width=1, space=c(1,1,2,1,1)) 
SET_NMS <- names(COMPILE[[1]])
SETS <- length(SET_NMS)
GENE_NMS <- names(COMPILE[[1]][[1]])
GENES <- length(GENE_NMS)
########################################################
## Compile Results for Plotting ########################
########################################################

## Make Array w/ all P-Vals
P_COMP <- array(,dim=c(length(GENE_NMS),length(SET_NMS)))
colnames(P_COMP) <- SET_NMS
rownames(P_COMP) <- GENE_NMS
for (i in 1:length(SET_NMS)) { P_COMP[,i] <- COMPILE[[1]][[i]] }

## Make Array w/ all Num_Vars
NV_COMP <- array(,dim=c(length(GENE_NMS),length(SET_NMS)+2))
colnames(NV_COMP) <- c(SET_NMS[1:4], paste(c("C","C","R","R"),SET_NMS[5:6],sep="_")[c(1,3,2,4)] )
rownames(NV_COMP) <- GENE_NMS
for (i in 1:4) { NV_COMP[,i] <- COMPILE[[5]][[i]] }
NV_COMP[,5:6] <- COMPILE[[5]][[5]]
NV_COMP[,7:8] <- COMPILE[[5]][[6]]
NV_COMP[which(is.na(NV_COMP))] <- 0

## Make Array of Rho Values
RO_COMP <- array(,dim=c(length(GENE_NMS),length(SET_NMS)))
colnames(RO_COMP) <- SET_NMS
rownames(RO_COMP) <- GENE_NMS
for (i in 1:length(SET_NMS)) { RO_COMP[,i] <- COMPILE[[6]][[i]] }

## Where does actual test rank in permuted p-vals?
G_Rank <- array(,c(GENES,SETS))
colnames(G_Rank) <- SET_NMS
rownames(G_Rank) <- GENE_NMS
for (g in 1:GENES) {
	for (s in 1:SETS) {
		G_Rank[g,s] <- length(which(COMPILE$P_Permuted[[s]][g,]<COMPILE$`P_SKAT-O`[[s]][g]))
	}	
}

## Scatter of P-Values and Counts
jpeg(paste(PathToBaseDir,"/Scatter_Pvals.jpeg",sep=""),height=1600,width=2400,pointsize=32)
YLIM <- c(0,max(-log10(P_COMP),na.rm=T))
XLIM <- c(1,GENES)
plot(0,0,type="n",xlim=XLIM,ylim=YLIM, xaxt="n", xlab="Genes", ylab="-log10(p-val)", main="P-Values for SKAT Tests by Gene")
for (i in 1:SETS) { points(1:GENES,-log10(t(P_COMP[,i])), col=SET_COLS[i], pch=20, type="b") }
legend(1,max(YLIM), legend=SET_NMS, col=SET_COLS, pch=20, title="Collapsed Set")
axis(1,at=1:GENES,labels=GENE_NMS, tick=T, las=2)
dev.off()

## Barplot of P-Values and Counts
jpeg(paste(PathToBaseDir,"/Barplot_Pvals_Counts.jpeg",sep=""),height=1600,width=2400,pointsize=32)
par(mfrow=c(2,1))
barplot(t(-log10(P_COMP)),beside=T,col=SET_COLS,legend=T,width=1,space=c(0,1), xlim=c(0,(SETS+1)*GENES) )
barplot(t(NV_COMP[,1:4]),beside=T,col=SET_COLS[1:4],legend=F,width=1,space=c(0,3), xlim=2+c(0,(SETS+1)*GENES) )
barplot(t(NV_COMP[,6:5]),beside=F,col=c("purple1","purple3"),width=1,space=c(7,rep(6,GENES-1)), add=T, axisnames=F )
barplot(t(NV_COMP[,8:7]),beside=F,col=c("gold1","gold3"),width=1,space=c(8,rep(6,GENES-1)), add=T, axisnames=F )
dev.off()








########################################################
## End of Doc ##########################################
########################################################