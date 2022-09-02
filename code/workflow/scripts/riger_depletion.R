libA_norm=read.table("normcounts_humanA_Mel1and2.txt",header = T,stringsAsFactors=F) 
libB_norm=read.table("normcounts_humanB_Mel1and2.txt",header = T,stringsAsFactors=F) 

# Give columns the same names in A and B libraries
temp = libA_norm
libName ="humanA_"
indsToReplace = grep(libName,colnames(temp))
for (i in 1:length(indsToReplace)) {
	colnames(temp)[indsToReplace[i]] <- sub(libName,"",colnames(temp)[indsToReplace[i]])
}
tempA = temp

temp = libB_norm
libName ="humanB_"
indsToReplace = grep(libName,colnames(temp))
for (i in 1:length(indsToReplace)) {
	colnames(temp)[indsToReplace[i]] <- sub(libName,"",colnames(temp)[indsToReplace[i]])
}
tempB = temp

# Unite A and B library
all_norm = rbind(tempA,tempB)

tdiff_d11_r1 = all_norm$Day11_T_12hr_Mel1 - all_norm$Day11_NoT_Mel1
tdiff_d17_r1 = all_norm$Day17_T_6hr_Mel1 - all_norm$Day17_NoT_Mel1
tdiff_d17_6hr_r2 = all_norm$Day17_T_6hr_Mel2 - all_norm$Day17_NoT_Mel2
tdiff_d17_12hr_r2 = all_norm$Day17_T_12hr_Mel2 - all_norm$Day17_NoT_Mel2

depletionDiff_d17_minus_d6 = all_norm$Day6_NoT_Mel1 - all_norm$Day17_NoT_Mel1

# Neville's weighted sum RIGER
dataForRIGER = depletionDiff_d17_minus_d6
# Rank scores for RIGER. For ties, use max ranking.
rigerDF = data.frame(UID=all_norm$UID, gene_id=all_norm$gene_id, rawDiff=dataForRIGER, rankedList=rank(-dataForRIGER, ties.method = "max"))
write.table(rigerDF,file="rigerDF_for_depletion.txt",row.names=F,quote=FALSE,sep="\t")
maxRankPercentage = 1
maxRank = round(maxRankPercentage/100*nrow(rigerDF))

# Calculates the RIGER Weighted Sum metric for a rigerDF containing 1) gene_id and 2) rankedList
# Make sure the rankedList is calculated using the R function "rank" (and not sort/order)
CalcWeightedSum <- function(rigerDF,maxRank=length(unique(rigerDF$gene_id))) {
	rankedDF=rigerDF[order(rigerDF$rankedList,decreasing = F),]
	
	# Eliminate sgRNAs where no score for ranking exists
	x=which(is.na(rankedDF$rawDiff))
	rankedDF=rankedDF[-x,]
	
	geneList <- NULL
	rigerScoreList <- NULL
	genesToConsider = unique(rankedDF$gene_id[1:maxRank])
	for (i in 1:length(genesToConsider)) {
		curGene = genesToConsider[i]
		top2Rank = sort(subset(rankedDF,gene_id==curGene)$rankedList)[1:2]
		geneList = c(geneList,as.character(curGene))
		# RIGER WeightedSum computation from: 
		# http://www.broadinstitute.org/cancer/software/GENE-E/extensions.html
		rigerScoreList= c(rigerScoreList, (0.25*top2Rank[1] + 0.75*top2Rank[2]))
	}
	outDF <- data.frame(gene_id=geneList, riger=rigerScoreList)
	outDF <- outDF[order(outDF$riger, decreasing=F),]
	return(outDF)
}		
		

RigerMC <- function(rigerDF, maxRank, numMC) {
	rigerValList=list()

	for (i in 1:numMC) {
		# Randomly permute ranks
		rigerDF$rankedList = sample.int(nrow(rigerDF))

		# Get RIGER stats
		curDF = CalcWeightedSum(rigerDF=rigerDF,maxRank=maxRank)
		rigerValList[[i]]= curDF$riger
	}
	return(rigerValList)
}

# adjustedMaxRank is the maximal rank of the best scoring guide to still have a RIGER rank < maxRank
# adjustedMaxRank = maxRank / 0.25
outWS <- CalcWeightedSum(rigerDF,maxRank=1200)
rvl10_1perc = RigerMC(rigerDF,maxRank=1200,numMC=10)
gene200 <- NULL; for (i in 1:10) { gene200[i]=rvl10_1perc[[i]][200];}

# Calculate 1% FDR for 20K genes --> top 200 genes
# Compute the RIGER score for #200 gene in each list
# Everything above that is >1% FDR for RIGER WS rank

# Calculate empirical CDF on all MC values computed above
ecdfFunc = ecdf(unlist(rvl10_1perc))
outWS$sigFDR = ecdfFunc(outWS$riger)

# Calculate significance based on random draws from [1,numSgs]
numSg = length(dataForRIGER)
drawOne=round(runif(1000000,min=0,max=numSg))
drawTwo=round(runif(1000000,min=0,max=numSg))
indsWhereOneIsLowerRank = which(drawOne>drawTwo)
x=drawOne[indsWhereOneIsLowerRank]*0.75 + drawTwo[indsWhereOneIsLowerRank]*0.25

ecdfFunc2 = ecdf(x)
outWS$sigRD = ecdfFunc2(outWS$riger)
outWS$RIGER_rank <- 1:nrow(outWS)
write.table(outWS,file="rigerOUTPUT_WS_mcSig10_depletion.txt",row.names=F,quote=FALSE,sep="\t")

# Compare to average rank of ribosomal genes
tWS = read.table("rigerOUTPUT_WS_mcSig10.txt", header = T)
tWS$RIGER_rank <- 1:nrow(tWS)
allRibGeneNames=unique(libA_norm$gene_id[grep("^RP[S|L]",libA_norm$gene_id)])

dataForRIGER = tdiff_d11_r1
# Rank scores for RIGER. For ties, use max ranking.
rigerDF = data.frame(UID=all_norm$UID, gene_id=all_norm$gene_id, rawDiff=dataForRIGER, rankedList=rank(-dataForRIGER, ties.method = "max"))
maxRankPercentage = 100
maxRank = round(maxRankPercentage/100*nrow(rigerDF))
x=is.na(rigerDF$rawDiff)
tWS2 <- CalcWeightedSum(rigerDF,maxRank=length(unique(rigerDF[!x,"gene_id"])))
tWS2$RIGER_rank <- 1:nrow(tWS2)
write.table(tWS2,file="rigerOUTPUT_tWS2_all20Kranks.txt",row.names=F,quote=FALSE,sep="\t")


y=subset(tWS2, gene_id %in% allRibGeneNames);

tWS3 = tWS2
tWS3$fracMut <- NA
for (i in 1:1000) {
	gene = tWS3$gene_id[i]
	x=getProfileData(cbiop,gene,variable,cases)
	tWS3$fracMut[i] = sum(!is.na(x[,1]))/length(x[,1])
}


# Below doesn't work
cosmicFilenames = list.files("/Users/neville/Downloads/Cosmic_data/")

for (i in 2:length(cosmicFilenames)) {
	x = read.table(paste0("/Users/neville/Downloads/Cosmic_data/",cosmicFilenames[i]),sep = '\t')
	tabNew = as.data.frame(table(x$V1))
	q=merge(tabOld,tabNew,by='Var1')
	tabNew=data.frame(Var1=q[,1], Freq=q[,2]+q[,3])
	tabOld=tabNew
}

# Gene names only
x=read.table("/Users/neville/Downloads/CosmicGeneNameOnly.tsv",header=T, colClasses = c("character"), sep='\t')
allDash=grep("_",x$Gene.name)
q=strsplit(x$Gene.name[allDash],split='_')
q2=unlist(lapply(q, `[[`, 1))

x2=x
x2$Gene.name[allDash] <- q2
tab=as.data.frame(table(x2))
indOrder = order(tab[,2], decreasing=T)
tab = tab[indOrder,]
write.table(tab,file=paste0('Cosmic_v80_AllMutationCounts',Sys.Date(),'.txt',sep='_'),row.names=F,quote=FALSE,sep="\t")

tWS4 = merge(tWS2,tab,by.x='gene_id',by.y='Gene')


geneCensus = read.csv('Census_allSat Apr 29 21-37-00 2017.csv', stringsAsFactors=F)
tsg=subset(geneCensus, Role.in.Cancer=='TSG')
tWS5 = merge(tWS2,tsg,by.x='gene_id',by.y='Gene.Symbol')
