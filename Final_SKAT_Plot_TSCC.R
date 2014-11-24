## Compile, Analyze, Plot SKATout.Rdata for all Runs ##
## August 27, 2014 ##
## Kristopher Standish ##

#############################################
## LOAD DATA ################################
#############################################

# Set Date
DATE <- "20140910"

# Set up Path to Data
PathToData <- paste("/projects/janssen/RV_Rare_Variant/",DATE,"_SKAT_Results/",sep="")
PathToPlot <- paste("/projects/janssen/RV_Rare_Variant/",DATE,"_SKAT_Results/",sep="")

# Get list of all files in Directory
File_List <- list.files(PathToData)
print(paste(length(File_List),"Rdata Files to Load"))

ALL <- list() ; start <- proc.time()
Loop <- 0
for (file in File_List) {
	Loop <- Loop + 1
	Load_File <- paste(PathToData,file,sep="")
	File_Name <- strsplit(file,"_LT8",fixed=T)[[1]][1]
	load(Load_File)
	ALL[[File_Name]] <- COMPILE
	if ( Loop%%10 == 0 ) {
		print(paste( "Loaded",Loop,"-",round((proc.time()-start)[3],1) ))
	}
}

#############################################
## COMPILE DATA INTO TABLE ##################
#############################################

## Pull out Info about which Sets were tested
P_COLS <- c("chartreuse1","dodgerblue1","firebrick1")
NV_COLS <- c("chartreuse1","gold1","firebrick1")
SET_NMS <- names(COMPILE[[1]])[c(1,2,4)]
SETS <- length(SET_NMS)

## Create variables in which to compile results
ALL_GENE_LIST <- c()
ALL_P_COMP <- array( ,dim=c(0,SETS) )
colnames(ALL_P_COMP) <- SET_NMS
ALL_NV_COMP <- array( ,dim=c(0,6) )
colnames(ALL_NV_COMP) <- c(paste("Full",c("C","R"),sep="_"),paste("Damg",c("C","R"),sep="_"),paste("Func",c("C","R"),sep="_"))
ALL_RO_COMP <- array( ,dim=c(0,SETS) )
colnames(ALL_RO_COMP) <- SET_NMS

## Loop through all different SKAT runs
for (rr in 1:length(ALL)) {
# for (rr in 1:5) {
	run_name <- names(ALL)[rr]
	COMPILE <- ALL[[rr]]

	## Pull out Gene names for this Run
	GENE_NMS <- names(COMPILE[[1]][[1]])
	GENES <- length(GENE_NMS)

	## Make Array w/ all P-Vals
	P_COMP <- array(,dim=c(length(GENE_NMS),length(SET_NMS)))
	colnames(P_COMP) <- SET_NMS
	rownames(P_COMP) <- GENE_NMS
	for (i in 1:length(SET_NMS)) { P_COMP[,i] <- COMPILE[["P_SKAT-O"]][[SET_NMS[i]]] }
	if ( all(is.na(P_COMP[,"Full"])) ) { print(paste("ALL MISSING P-VALUES FOR",run_name)) }

	## Make Array w/ all Num_Vars
	NV_COMP <- array(,dim=c(length(GENE_NMS),6))
	# colnames(NV_COMP) <- c(SET_NMS[1],paste(SET_NMS[2],c("C","R"),sep="_")) # c(SET_NMS[1:4], paste(c("C","C","R","R"),SET_NMS[5:6],sep="_")[c(1,3,2,4)] )
	colnames(NV_COMP) <- c(paste("Full",c("C","R"),sep="_"),paste("Damg",c("C","R"),sep="_"),paste("Func",c("C","R"),sep="_")) # c(SET_NMS[1:4], paste(c("C","C","R","R"),SET_NMS[5:6],sep="_")[c(1,3,2,4)] )
	rownames(NV_COMP) <- GENE_NMS
	NV_COMP[,grep("Full",colnames(NV_COMP))] <- COMPILE$Num_Vars$Full_All
	NV_COMP[,grep("Damg",colnames(NV_COMP))] <- COMPILE$Num_Vars$Damg_All
	NV_COMP[,grep("Func",colnames(NV_COMP))] <- COMPILE$Var_Counts[,grep("Func",colnames(COMPILE$Var_Counts))]
	NV_COMP[which(is.na(NV_COMP))] <- 0

	## Make Array of Rho Values
	RO_COMP <- array(,dim=c(length(GENE_NMS),length(SET_NMS)))
	colnames(RO_COMP) <- SET_NMS
	rownames(RO_COMP) <- GENE_NMS
	for (i in 1:length(SET_NMS)) { RO_COMP[,i] <- COMPILE[["Rho_Val"]][[SET_NMS[i]]] }

	## Compile results from separate runs into 1 table
	ALL_GENE_LIST <- c(ALL_GENE_LIST, GENE_NMS)
	ALL_P_COMP <- rbind(ALL_P_COMP, P_COMP)
	ALL_NV_COMP <- rbind(ALL_NV_COMP, NV_COMP)
	ALL_RO_COMP <- rbind(ALL_RO_COMP, RO_COMP)

} ## Close Run loop (rr)

## Write some of the tables so I don't have to compile this shiz again
write.table(ALL_GENE_LIST,paste(PathToPlot,"ALL_Gene_List.txt",sep=""), sep="\t",row.names=F,col.names=F,quote=F)
write.table(ALL_P_COMP,paste(PathToPlot,"ALL_P_Comp.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F)
write.table(ALL_NV_COMP,paste(PathToPlot,"ALL_NV_Comp.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F)
write.table(ALL_RO_COMP,paste(PathToPlot,"ALL_RO_Comp.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F)
# ALL_GENE_LIST <- read.table(paste(PathToPlot,"ALL_Gene_List.txt",sep=""), sep="\t",header=T)
# ALL_P_COMP <- read.table(paste(PathToPlot,"ALL_P_Comp.txt",sep=""), sep="\t",header=T)
# ALL_NV_COMP <- read.table(paste(PathToPlot,"ALL_NV_Comp.txt",sep=""), sep="\t",header=T)
# ALL_RO_COMP <- read.table(paste(PathToPlot,"ALL_RO_Comp.txt",sep=""), sep="\t",header=T)

#############################################
## PLOT THIS SHIZ ###########################
#############################################

## Make QQ Plot for this Shiz
print("Making QQ-Plot")
X_MAX <- ceiling( -log10(1/nrow(ALL_P_COMP)) )
Y_MAX <- ceiling( -log10(min(ALL_P_COMP,na.rm=T)) )
MAX <- max(X_MAX,Y_MAX)

png(paste(PathToPlot,"ALL_",DATE,"_1-QQ_Plot.png",sep=""), height=1600, width=1600, pointsize=30)
plot(0,0,type="n",xlim=c(0,MAX),ylim=c(0,MAX),xlab="Expected",ylab="Observed",main="Q-Q for SKAT Tests on Collapsed Sets")
abline(0,1,lty=2,lwd=4,col="grey10")
abline(h=1:8, lty=2, lwd=1, col="grey50")
abline(v=1:8, lty=2, lwd=1, col="grey50")
for (i in 1:ncol(ALL_P_COMP)) {
  if ( length(which( !is.na(ALL_P_COMP[,i]) )) > 0 ) {
    OBS <- sort(ALL_P_COMP[,i],na.last=NA)
    print(paste("length(OBS)=",length(OBS)))
    # EXP <- 1:N_TESTS/N_TESTS
    EXP <- 1:length(OBS)/length(OBS)
    print(paste("length(EXP)=",length(EXP)))
    points(-log10(EXP),-log10(OBS),col=P_COLS[i],pch="+")
    text(.5,4-.2*i, labels=paste("# Genes:",length(OBS)), col=P_COLS[i])
  }
}
legend(0.7*MAX,0.3*MAX,legend=SET_NMS,pch="+",col=P_COLS)
dev.off()

## Scatter P-Vals for "Full" vs "Damg" sets
print("Full vs Damg vs Full_All")
X_MAX <- ceiling( -log10(1/nrow(ALL_P_COMP)) )
Y_MAX <- ceiling( -log10(min(ALL_P_COMP,na.rm=T)) )
MAX <- max(X_MAX,Y_MAX)

png(paste(PathToPlot,"ALL_",DATE,"_1b-Full_v_Damg.png",sep=""), height=2000, width=2000, pointsize=30)
par(mfrow=c(2,2))
# Full_ALL v Full_R
plot(0,0,type="n",xlim=c(0,MAX),ylim=c(0,MAX),xlab="Full (All)",ylab="Full (Rare)",main="P-Values: Full (All) vs Full (Rare) Sets")
abline(0,1,lty=2,lwd=2,col="grey50")
WHICH <- which(!is.na(ALL_P_COMP[,"Full_All"]))
points(-log10(ALL_P_COMP[WHICH,"Full_All"]),-log10(ALL_P_COMP[WHICH,"Full"]), pch="+", col="cadetblue1")
# Empty
plot(0,0,type="n",xlim=c(0,MAX),ylim=c(0,MAX),xlab="",ylab="")
# Full_ALL v Damg
plot(0,0,type="n",xlim=c(0,MAX),ylim=c(0,MAX),xlab="Full (All)",ylab="Damg (All)",main="P-Values: Full (All) vs Damg Sets")
abline(0,1,lty=2,lwd=2,col="grey50")
WHICH <- which(!is.na(ALL_P_COMP[,"Damg_All"]))
points(-log10(ALL_P_COMP[WHICH,"Full_All"]),-log10(ALL_P_COMP[WHICH,"Damg_All"]), pch="+", col="purple2")
# Full_R v Damg
plot(0,0,type="n",xlim=c(0,MAX),ylim=c(0,MAX),xlab="Full (Rare)",ylab="Damg (All)",main="P-Values: Full (Rare) vs Damg Sets")
abline(0,1,lty=2,lwd=2,col="grey50")
WHICH <- which(!is.na(ALL_P_COMP[,"Damg_All"]))
points(-log10(ALL_P_COMP[WHICH,"Full"]),-log10(ALL_P_COMP[WHICH,"Damg_All"]), pch="+", col="gold1")
dev.off()

## Make Scatter of P-Values for Top Hits
n_FULL <- length(which( !is.na(ALL_P_COMP[,"Full"]) ))
n_FULL_ALL <- length(which( !is.na(ALL_P_COMP[,"Full_All"]) ))
n_DAMG <- length(which( !is.na(ALL_P_COMP[,"Damg_All"]) ))
THRSH <- 100/nrow(ALL_P_COMP)
BEST_P_FULL <- ALL_P_COMP[ which(ALL_P_COMP[,"Full"] < quantile(ALL_P_COMP[,"Full"],THRSH,na.rm=T) ), ]
BEST_P_FULL_ALL <- ALL_P_COMP[ which(ALL_P_COMP[,"Full_All"] < quantile(ALL_P_COMP[,"Full_All"],THRSH,na.rm=T) ), ]
BEST_P_DAMG <- ALL_P_COMP[ which(ALL_P_COMP[,"Damg_All"] < quantile(ALL_P_COMP[,"Damg_All"],THRSH,na.rm=T) ), ]
BEST_P_NAMES <- sort( Reduce(union,list(rownames(BEST_P_FULL),rownames(BEST_P_DAMG),rownames(BEST_P_FULL_ALL))) )
BEST_P <- ALL_P_COMP[BEST_P_NAMES,]
BEST_NV <- ALL_NV_COMP[BEST_P_NAMES,]
BEST_RO <- ALL_RO_COMP[BEST_P_NAMES,]
 # Scatter of P-Values
print("Plotting P-Values")
PLOT_WD <- 800+25*nrow(BEST_P)
png(paste(PathToPlot,"ALL_",DATE,"_2-P_Scatter.png",sep=""), height=1600, width=PLOT_WD, pointsize=30)
YLIM <- c(0,max(-log10(BEST_P),-log10(.05/n_FULL),na.rm=T))
XLIM <- c(1,nrow(BEST_P))
plot(0,0,type="n",xlim=XLIM,ylim=YLIM, xaxt="n", xlab="",ylab="-log10(p-val)", main="P-Values for SKAT Tests by Gene")
legend(1,max(YLIM), legend=SET_NMS, col=P_COLS, pch="+", title="Collapsed Set")
axis(1,at=1:nrow(BEST_P),labels=rownames(BEST_P), tick=T, las=2)
abline(h=-log10(.05/n_FULL), col=P_COLS[1], lty=2, lwd=3) ; text(0.95*XLIM[2],.05-log10(.05/n_FULL), labels="Bonf - Full", col=P_COLS[1])
abline(h=-log10(.05/n_FULL_ALL), col=P_COLS[2], lty=2, lwd=3) ; text(0.95*XLIM[2],.05-log10(.05/n_FULL_ALL), labels="Bonf - Full_All", col=P_COLS[2])
abline(h=-log10(.05/n_DAMG), col=P_COLS[3], lty=2, lwd=3) ; text(0.95*XLIM[2],.05-log10(.05/n_DAMG), labels="Bonf - Damg", col=P_COLS[3])
for (i in 1:SETS) { points(1:nrow(BEST_P),-log10(t(BEST_P[,i])), col=P_COLS[i], pch="+", type="b") }
dev.off()

## Barplot of Counts
print("Plotting # Variants Used in Tests")
PLOT_WD <- 800+25*nrow(BEST_P)
png(paste(PathToPlot,"ALL_",DATE,"_3-NV_Barplot.png",sep=""), height=1600, width=PLOT_WD, pointsize=30)
par(mfrow=c(3,1))
barplot(t(BEST_NV[,grep("Full",colnames(BEST_NV))]),beside=F,col=c(NV_COLS[1],gsub("1","4",NV_COLS[1])),legend=F,width=.8,space=.25, ylim=c(0,1.2*max(BEST_NV[,grep("Full",colnames(BEST_NV))],na.rm=T)), las=2, main="Variants Included in Full Test", ylab="# Variants")
legend(0,1.4*max(rowSums(BEST_NV[,grep("Full",colnames(BEST_NV))]),na.rm=T),legend=c("Rare","Common"),fill=c(gsub("1","4",NV_COLS[1]),NV_COLS[1]))
abline(h=8000, lty=2, col="grey50")
barplot(t(BEST_NV[,grep("Func",colnames(BEST_NV))]),beside=F,col=c(NV_COLS[2],gsub("1","4",NV_COLS[2])),legend=F,width=.8,space=.25, ylim=c(0,1.2*max(BEST_NV[,grep("Func",colnames(BEST_NV))],na.rm=T)), las=2, main="Variants in Exonic/UTR Regions", ylab="# Variants")
legend(0,1.4*max(rowSums(BEST_NV[,grep("Func",colnames(BEST_NV))]),na.rm=T),legend=c("Rare","Common"),fill=c(gsub("1","4",NV_COLS[2]),NV_COLS[2]))
barplot(t(BEST_NV[,grep("Damg",colnames(BEST_NV))]),beside=F,col=c(NV_COLS[3],gsub("1","4",NV_COLS[3])),legend=F,width=.8,space=.25, ylim=c(0,1.2*max(BEST_NV[,grep("Damg",colnames(BEST_NV))],na.rm=T)), las=2, main="Nonsynonymous Variants Included in Test ('Damg' Set)", ylab="# Variants")
legend(0,1.4*max(rowSums(BEST_NV[,grep("Damg",colnames(BEST_NV))]),na.rm=T),legend=c("Rare","Common"),fill=c(gsub("1","4",NV_COLS[3]),NV_COLS[3]))
dev.off()

## Scatter P-Vals vs N_Vars (both Sets)
png(paste(PathToPlot,"ALL_",DATE,"_14-Pval_v_Count.png",sep=""), height=800, width=2400, pointsize=30)
par(mfrow=c(1,3))
# plot(0,0,type="n",xlab="# Variants",ylab="P-Values (-log10)",main="P-Values vs Number of Rare Variants in Gene")
# plot(log10(ALL_NV_COMP[,"Full"]),-log10(ALL_P_COMP[,"Full"]), pch="+", col=P_COLS[1],xlab="# Variants (log10)",ylab="P-Values (-log10)",main="P-Values vs Number of Rare Variants in Gene")
# plot(log10(ALL_NV_COMP[,"Full_All"]),-log10(ALL_P_COMP[,"Full_All"]), pch="+", col=P_COLS[2],xlab="# Variants (log10)",ylab="P-Values (-log10)",main="P-Values vs Total Number of Variants in Gene")
plot(log10(ALL_NV_COMP[,"Full_R"]),-log10(ALL_P_COMP[,"Full"]), pch="+", col=P_COLS[1],xlab="# Variants (log10)",ylab="P-Values (-log10)",main="P-Values vs Number of Rare Variants in Gene")
plot(log10(rowSums(ALL_NV_COMP[,c("Full_C","Full_R")])),-log10(ALL_P_COMP[,"Full_All"]), pch="+", col=P_COLS[2],xlab="# Variants (log10)",ylab="P-Values (-log10)",main="P-Values vs Total Number of Variants in Gene")
WHICH <- which(!is.na(ALL_P_COMP[,"Damg_All"]))
# plot(0,0,type="n",xlab="# Variants",ylab="P-Values (-log10)",main="P-Values vs Total Nonsynonymous Variants in Gene")
plot(log10(rowSums(ALL_NV_COMP[WHICH,c("Damg_C","Damg_R")])),-log10(ALL_P_COMP[WHICH,"Damg_All"]), pch="+", col=P_COLS[3],xlab="# Variants (log10)",ylab="P-Values (-log10)",main="P-Values vs Nonsynonymous Variants in Gene")
dev.off()


#############################################
## END OF DOC ###############################
#############################################
