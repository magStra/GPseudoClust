#Computing cluster allocations using alternative algorithms
#and comparing to true cluster allocations
#First run sumClustSim.R to find summary clusterings for results obtained by GPseuodClust
#mclust
simData <- read.csv("simDataClust.csv",header=F)
options(stringsAsFactors = F)
library(mclust)
library(cluster)

simGeneNames <- 1:dim(simData)[1]

clusters_mclust <- Mclust(simData)$classification

#k-means
library(fpc)
clusters_KMeans <- pamk(simData,1:15)$pamobject$clustering

clusters_H <- hclust(dist(simData))
library(clValid)
a <- clValid(simData,2:10,"hierarchical",validation="internal")
summary(a)
clusters_H <- cutree(clusters_H,5)


library(SLICER)
library(lle)
traj = t(simData)
genes <- select_genes(traj)
traj <- traj[,genes]
k = select_k(traj, kmin=5,kmax=50)
traj_lle = lle(traj, m=2, k)$Y
traj_graph = conn_knn_graph(traj_lle,5)
ends = find_extreme_cells(traj_graph, traj_lle)
start = ends[2]
cells_orderedSLICER = cell_order(traj_graph, start)
pos <- 1:60
pseudoTimesSLICER <- sort(process_distance(traj_graph, start))
assign_branches(traj_graph,start=ends[2])
###SLICER recognises there is no branching
simSLICER <- simData[,cells_orderedSLICER]
write.table(simSLICER,file = "simSLICER.csv",sep=",",row.names = F,col.names = F)
write.table(pseudoTimesSLICER,file = "pseudoTimesSimSLICER.csv",sep=",",row.names = F,col.names = F)

#after ordering, clustering performed by GPClust method, see GPClustSimulations.ipynb and save as "GPClustSim1Slicer.csv"
GP_Clust <- unlist(read.table("GPClustSim1Slicer.csv",sep=",",header=F))


# ordering and clustering with monocle2
#Uncomment to run, alternatively (as there are namespace problems with some
#set-ups and versions, we provide the cluster allocations: )

library(monocle)
pd <- data.frame(1:60)
pd.pd <- colnames(simData)
rownames(pd)=colnames(simData)
pd <- as(pd, "AnnotatedDataFrame")

fd <- data.frame(1:52)
fd$gene_short_name <- sapply(1:52,toString)
fd <- fd[,2,drop=F]
rownames(fd) <- rownames(simData)
fd <- as(fd, "AnnotatedDataFrame")
cds <- newCellDataSet(as.matrix(simData), phenoData = pd, featureData = fd,expressionFamily = uninormal())
cds <-  reduceDimension(cds,pseudo_expr = 0,verbose = T,reduction_method="ICA",norm_method = "none")
cds <- orderCells(cds,num_paths = 1)
full_model_fits <- fitModel(cds)
expression_curve_matrix <- responseMatrix(full_model_fits)#' This function fits a vector generalized additive model (VGAM) from the VGAM package for each gene in a CellDataSet.
#' By default, expression levels are modeled as smooth functions of the Pseudotime value of each
#' cell.

distMonocle <- as.dist(1-cor(t(expression_curve_matrix))/2)
sil <- c()
clusters_Monocle<-list()
for (j in 2:15){
  clusters_Monocle[[j-1]] <-  clusterGenes(expression_curve_matrix,k=j)$cluster
  sil = c(sil,mean(silhouette(clusters_Monocle[[j-1]],distMonocle)[,3]))}
clusters_Monocle = clusters_Monocle[[which.max(sil)]]

#deLorean and GPClust
#uncomment to run
# library(DeLorean)
# pd <- data.frame(cell =as.factor(sapply(1:60,toString)),obstime= c(rep(1,20),rep(2,20),rep(3,20)),
#                 capture = as.factor(c(rep(1,20),rep(2,20),rep(3,20))) )
# fd <- data.frame(gene = as.factor(simGeneNames))
# sim.expr <- as.matrix(simData)
# rownames(fd) <- simGeneNames
# rownames(pd) <- sapply(1:60,toString)
# rownames(sim.expr) <- simGeneNames
# colnames(sim.expr) <- sapply(1:60,toString)
# dl <- de.lorean(sim.expr, fd, pd)
# dl <- estimate.hyper(
#   dl,
#   sigma.tau=1,
#   length.scale=3,
#   model.name='lowrank-sizes')
# SimDL <- fit.dl(dl, method='vb')
# SimDL <- examine.convergence(SimDL)
# plot(SimDL, type='Rhat')
# sample.iter=SimDL$best.sample
# SimDL$samples.l$tau
# pseudoTimesDL <- SimDL$samples.l$tau$tau[SimDL$samples.l$tau$iter == sample.iter]
# pseudoTimeOrder <- sort(pseudoTimesDL,index.return=T)
# write.table(pseudoTimeOrder$x,file="deLoreanTSim.csv",sep=",",row.names=F,col.names=F)
# write.table(simData[,pseudoTimeOrder$ix],file="deLoreanSim.csv",sep=",",row.names=F,col.names=F)

#ordering performed by GPClust method, see GPClustSimulations.ipynb and save as "GPClustSim1.csv"
GP_Clust1 <- unlist(read.table("GPClustSim1DL.csv",sep=",",header=F))

#import clusters found by SIMLR
clustersSIMLR <- unlist(read.table("../postProcessing/sim1_SIMLR.csv"))
#load GPseudoClust summary clusterings
GPseudoClust_lmkk <- read.table('simlmkk.csv')
load("npSumClustSim.rda")

library(ClusterR)
true_clusters <- c(rep(1,8),rep(2,4),rep(3,12),rep(4,16),rep(5,12))

externalValidation <- function(clusterAllocations)
{
  output <- rep(0,3)
  output[1] <- adjustedRandIndex(clusterAllocations,true_clusters)
  output[2] <- external_validation(true_clusters, clusterAllocations, method = "fowlkes_mallows_index",
                                   summary_stats = F)
  output[3] <- external_validation(true_clusters, clusterAllocations, method = "nmi",
                                   summary_stats = F)
  output
}


allAllocs <- cbind(clusters_mclust,clustersSIMLR,clusters_KMeans,clusters_H, GP_Clust,GP_Clust1,clusters_Monocle,
                   GPseudoClust_lmkk,GPseudoClust_PY,GPseudoClust_DP,GPseudoClust_mean)

evalMatSim <- apply(allAllocs,2,externalValidation)
colnames(evalMatSim) <- c('mcl','SIMLR','PAM','hier.','SL+GCl','De+GCl','Mon 2','GPs+lmkk','GPs+PY','GPs+DP','GPs+mean')
write.table(evalMatSim,file="clusterCompSim1.csv",sep=",",row.names=F)
