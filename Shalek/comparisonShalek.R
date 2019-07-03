options(stringsAsFactors = F)
#read preprocessed Shalek data set
ShalekData <- read.table("../../preprocessedScRNASeq/Shalek.csv",sep=",",header=T)
ShalekGeneNames <- ShalekData[,1]
ShalekData <- ShalekData[,126:308]#excluding gene names, and capture times 0h, 1h
library(mclust)
library(cluster)
# 
clusters_mclust <- Mclust(ShalekData,G=4:20)$classification

library(fpc)
clusters_KMeans <- pamk(ShalekData,4:15)$pamobject$clustering

clusters_H <- hclust(dist(ShalekData))
plot(clusters_H)
library(clValid)
clComp <- clValid(ShalekData, 4:15, clMethods = "hierarchical", validation =
          "internal", maxitems = 600, metric = "euclidean", method = "complete",
        neighbSize = 10, annotation = NULL, 
         dropEvidence=NULL, verbose=FALSE)
summary(clComp)
clusters_H <- cutree(clusters_H,k=4)#from silhouette width with clValid

library(SLICER)
library(lle)
traj = t(ShalekData)
genes <- select_genes(traj)
traj <- traj[,genes]
k = select_k(traj, kmin=5,kmax=50)
traj_lle = lle(traj, m=2, k)$Y
traj_graph = conn_knn_graph(traj_lle,5)
ends = find_extreme_cells(traj_graph, traj_lle)
start = ends[2]
cells_orderedSLICER = cell_order(traj_graph, start)
pseudoTimesSLICER <- sort(process_distance(traj_graph, start))
assign_branches(traj_graph,start=ends[2])
                ###SLICER recognises there is no branching
ShalekSLICER <- ShalekData[,cells_orderedSLICER]
write.table(ShalekSLICER,file = "../scratch/ShalekSLICER.csv",sep=",",
            col.names = F,row.names=F)
write.table(pseudoTimesSLICER,file = "../scratch/pseudoTimesSLICER.csv",sep=",", col.names = F,row.names=F)
write.table(shQuote(ShalekGeneNames),file = "../scratch/ShalekGeneNames.csv",sep=",", col.names = F,row.names=F)

# ###clusters from ordering with SLICER and subsequent clustering with GPclust
# #ordering and GPclust, output of GPclust is sensitive to initialisation and may vary 
# a bit and need several restarts. The following is copied from 
GP_Clust <- c(3.,  2. , 4. , 2. , 3. , 1. , 1.  ,5. , 4. , 1. , 4.,  1. , 1. , 2.,  5. , 1.,  1.,  1.,
              2.,  1.,  1.,  4. , 2.,  5. , 2. , 3. , 3.,  3. , 3.,  2. , 1.,  1. , 2. , 3. , 1.,  3.,
              3, 1,  1,  6,  2, 2,  1,  1,  1,  2,  1,  2,  4,  1,  3,  1,  4,  1,
              1,  6,  3,  2,  4,  5,  1,  5,  3,  1,  4,  1,  2,  1,  2,  1,  1,  3,
              3,  1) # copied from Python notebook
# 
# # ordering with DeLorean and clustering with GPClust
# 
# library(DeLorean)
# pd <- data.frame(cell =as.factor(sapply(1:183,toString)),obstime= c(rep(2,65),rep(4,60),rep(6,58)),
#                 capture = as.factor(c(rep(2,65),rep(4,60),rep(6,58))) )
# fd <- data.frame(gene = as.factor(ShalekGeneNames))
# Shalek.expr <- as.matrix(ShalekData)
# rownames(fd) <- ShalekGeneNames
# rownames(pd) <- sapply(1:183,toString)
# rownames(Shalek.expr) <- ShalekGeneNames
# colnames(Shalek.expr) <- sapply(1:183,toString)
# dl <- de.lorean(Shalek.expr, fd, pd)
# dl <- estimate.hyper(
#   dl,
#   sigma.tau=1,
#   length.scale=3,
#   model.name='exact')
# ShalekDL <- fit.dl(dl, method='vb')
# ShalekDL <- examine.convergence(ShalekDL)
# plot(ShalekDL, type='Rhat')
# sample.iter=ShalekDL$best.sample
# ShalekDL$samples.l$tau
# pseudoTimesDL <- ShalekDL$samples.l$tau$tau[ShalekDL$samples.l$tau$iter == sample.iter]
# pseudoTimeOrder <- sort(pseudoTimesDL,index.return=T)
# write.csv(pseudoTimeOrder$x,file="../scratch/deLoreanTShalek.csv")
# write.csv(ShalekData[,pseudoTimeOrder$ix],file="/scratch/deLoreanShalek.csv")

GP_Clust1 <- c(2,  1,  1,  1,  6,  1,  1,  2,  1,  1,  1,  1,  1,  1,  4,  1,  1,  1,
               1,  1,  1,  1,  1,  2,  1,  5,  3,  2,  5,  1,  1,  1,  1,  3,  1,  3,
               2,  1,  1,  2,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  6,  1,  1,  1,
               1,  4,  5,  1,  2,  2,  1,  4,  3,  1,  1,  1,  1,  1,  1,  1,  1,  2,
               3,  1) #again copied from python output
adjustedRandIndex(GP_Clust,GP_Clust1)
GPclusters1 <- list()
for (j in 1:max(GP_Clust1)){
  GPclusters1[[j]] <- ShalekGeneNames[which(GP_Clust1 == j)]}

# ordering and clustering with monocle2
library(monocle)
load("../data/ShalekDeLorean.rda")#load data without cell size correction, as required by monocle2
ShalekSubset <- shalek.A.expr[ShalekGeneNames,colnames(ShalekData)]
pd <- shalek.A.cell.meta[colnames(ShalekData),]
rownames(pd)=colnames(ShalekData)
names(shalek.A.gene.meta)[1] <- "gene_short_name"
pd <- as(pd, "AnnotatedDataFrame")

fd <- shalek.A.gene.meta[xx,]
rownames(fd) <- ShalekGeneNames
fd <- as(fd, "AnnotatedDataFrame")
cds <- newCellDataSet(as.matrix(ShalekData), phenoData = pd, featureData = fd,expressionFamily = gaussianff())
diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~Media")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)#none are excluded as expected
cds <-  reduceDimension(cds,pseudo_expr = 0,verbose = T,norm_method = "none")
cds <- orderCells(cds)
full_model_fits <- fitModel(cds)
expression_curve_matrix <- responseMatrix(full_model_fits)#' This function fits a vector generalized additive model (VGAM) from the VGAM package for each gene in a CellDataSet.
#' By default, expression levels are modeled as smooth functions of the Pseudotime value of each
#' cell.
clusters <- clusterGenes(expression_curve_matrix,k=4)
plot_clusters(cds, clusters)
distMonocle <- as.dist(1-cor(t(expression_curve_matrix))/2)
sil <- c()
clusters_Monocle<-list()
for (j in 4:15){#4 modules in the data, so at least 4 clusters
  clusters_Monocle[[j-3]] <-  clusterGenes(expression_curve_matrix,k=j)$cluster
  sil = c(sil,mean(silhouette(clusters_Monocle[[j-3]],distMonocle)[,3]))}
clusters_Monocle = clusters_Monocle[[which.max(sil)]]

#in case of namespace error, load pre-computed clusters for monocle

load("../scratch/clustersMonocle.rda")

#import SIMLR clusters (SIMLRComp.m)
options(stringsAsFactors = F)
allocsSIMLR = read.table("../postProcessing/Shalek_SIMLR.csv",sep=",",header=F)$V1

#compute ARI between the different inferred cluster structures

library(ClusterR)
pairwiseARI <- function(listOfAllocations)
{
  l = length(listOfAllocations)
  output = 0.5*diag(l)
  for (j in 1:(l-1)){
    for (k in (j+1):l)
    {
      output[j,k] = external_validation(listOfAllocations[[k]],listOfAllocations[[j]],method = "adjusted_rand_index")
    }
  }
  output = output+t(output)
  output
}

listOfAllocations = list(clusters_mclust,allocsSIMLR,clusters_KMeans,clusters_H,GP_Clust,GP_Clust1,clusters_Monocle)
saveMat1 <- pairwiseARI(listOfAllocations)
write.table(saveMat1,"Shalek_ARI.csv",sep=",",row.names=F,col.names=F)
rownames(saveMat1) <- c('mcl','SIMLR','PAM','hier.','SL+GCl','De+GCl','Mon 2')
colnames(saveMat1) <- c('mcl','SIMLR','PAM','hier.','SL+GCl','De+GCl','Mon 2')
#plot heatmap of ARIs
heatmap(saveMat1, Colv = NA, Rowv = NA, scale="column")
#adding draws from GPseudoClust posterior
chains <- sample(1:96,8,replace=T)
samples <- sample(2001:4000,8,replace=T)
allocs <- c()
for (j in 1:8){
  fileName <- sprintf("Shalek_Results_Chain%d.csv",chains[j])
  A <- read.table(fileName ,sep=",",header=F,skip=chains[j])
  allocs <- rbind(allocs,A[1,2:75])
}
save(allocs,file="Shalek_randomPostSamples.rda")
write.table(allocs,"Shalek_randomPostSamples.csv",sep=",",row.names=F,col.names = F)

load("Shalek_randomPostSamples.rda")
allocs <- as.matrix(allocs)
allocsListGPs <- list()
for (j in 1:8)
{
  allocsListGPs[[j]] <- allocs[j,]
}
listOfAllocationsExt <- append(listOfAllocations,allocsListGPs)
saveMat2 <- pairwiseARI(listOfAllocationsExt)
write.table(saveMat2,"Shalek_ARIExt.csv",sep=",",row.names=F,col.names=F)
