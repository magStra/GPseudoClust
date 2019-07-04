#The raw data are available on Mendeley Data (http://dx.doi.org/10.17632/g2md5gbhz7.1)
#Stumpf PS et al. (2017). Stem Cell Differentiation as a Non-Markov Stochastic Process. Cell Systems 5, 268 - 282.e7.

preProcessStumpfData = function(cellLine,pathToData)
{
  if (cellLine == "E"){
  subS = 1:48
  Name = "E"}
  else{
    subS = 49:96
  Name = "R"
}
T0 = read.delim(paste(pathToData,"single-cell-gene-expression-data-of-mesc-differentiation-time-course-toward-the-neuronal-lineage./2015-01-08\ 0H\ sample\ Ct\ values.txt",sep=""))

T0 = T0[,2:97]
T1 = read.delim(paste(pathToData,"single-cell-gene-expression-data-of-mesc-differentiation-time-course-toward-the-neuronal-lineage./2015-01-08\ 24H\ sample\ Ct\ values.txt",sep=""))
geneNames = T1$Genes
write.csv(geneNames,"StumpfGeneNames.csv")
T1 = T1[,2:97]
T2 = read.delim(paste(pathToData,"single-cell-gene-expression-data-of-mesc-differentiation-time-course-toward-the-neuronal-lineage./2015-01-08\ 48H\ sample\ Ct\ values.txt",sep=""))
T2 = T2[,2:97]

T3 = read.delim(paste(pathToData,"single-cell-gene-expression-data-of-mesc-differentiation-time-course-toward-the-neuronal-lineage./2015-01-08\ 72H\ sample\ Ct\ values.txt",sep=""))
T3 = T3[,2:97]

T4 = read.delim(paste(pathToData,"single-cell-gene-expression-data-of-mesc-differentiation-time-course-toward-the-neuronal-lineage./2015-01-08\ 96H\ sample\ Ct\ values.txt",sep=""))
T4 = T4[1:96,2:97]

T5 = read.delim(paste(pathToData,"single-cell-gene-expression-data-of-mesc-differentiation-time-course-toward-the-neuronal-lineage./2015-01-08\ 120H\ sample\ Ct\ values.txt",sep=""))
#selecting cells with less than 50% NA
T5 = T5[1:96,2:97]

T6 = read.delim(paste(pathToData,"single-cell-gene-expression-data-of-mesc-differentiation-time-course-toward-the-neuronal-lineage./2015-01-08\ 168H\ sample\ Ct\ values.txt",sep=""))

T6 = T6[1:96,2:97]


data = cbind(T0[,subS],T1[,subS],T2[,subS],T3[,subS],T4[,subS],T5[,subS],T6[,subS])
data[data>28] <- 28
captureTimes <- c()
for (j in 1:7)
{captureTimes <- c(captureTimes,rep(j,48))}

ind <- c(which(geneNames == "Actb"),which(geneNames == "Gapdh"))
#remove the cells where Actb or Gapdh > 15
ind_keep <- data[ind[1],]<=15
data <- data[,ind_keep]
captureTimes <- captureTimes[ind_keep]
ind_keep <- data[ind[2],]<=15

#normalisation using loading controls
cellNames <- colnames(data)
colnames(data) <- c()
rownames(data) <- c()
data <- as.matrix(data)
uCT <- unique(captureTimes)
for (j in 1:7)
{
  data[,captureTimes=uCT[j]] <- data[,captureTimes=uCT[j]] - median(data[ind,captureTimes=uCT[j]])
}


#Compute quartiles of ET
ET <- colSums(data)
quartiles <- quantile(ET)
L <- quartiles[4] - quartiles[2]
threshLow <- quartiles[2]-2*L
threshUp <- quartiles[4]+2*L
ind_keep <- ET > threshLow
data <- data[,ind_keep]
cellNames <- cellNames[ind_keep]
captureTimes <- captureTimes[ind_keep]
ET <- colSums(data)
ind_keep <- ET < threshUp
data <- data[,ind_keep]
cellNames <- cellNames[ind_keep]
captureTimes <- captureTimes[ind_keep]
max(data)
min(data)

#a + b*min(data) = 28
#a + b*max(data) = 0
b = -28/(max(data)-min(data))
a = -b*max(data)

data <- b*data+a

#finally remove  Actb and Gapdh
data <- data[-ind,]
geneNames <- geneNames[-ind]

write.table(geneNames,"StumpfGeneNames.csv",sep=",",
            col.names=F,row.names=F)
write.table(data,paste("Stumpf_",Name,".csv",sep=""),sep=",",
            col.names=F,row.names=F)
write.table(cellNames,paste("StumpfCellNames_",Name,".csv",sep=""),sep=",",
            col.names=F,row.names=F)
write.table(captureTimes,paste("StumpfCT_",Name,".csv",sep=""),sep=",",
            col.names=F,row.names=F)
#
}

preProcessStumpfData("E","~/Documents/GPseudoRankCode/")
preProcessStumpfData("R","~/Documents/GPseudoRankCode/")
