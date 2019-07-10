simData <- read.csv("simDataClust2.csv",header=F)

library(mclust)
library(nonparametricSummaryPSM)

simGeneNames <- 1:dim(simData)[1]

clusters_mclust <- Mclust(simData)$classification
#PAM
library(fpc)
clusters_KMeans <- pamk(simData,1:15)$pamobject$clustering

clusters_H <- hclust(dist(simData))
library(clValid)
a <- clValid(simData,2:10,"hierarchical",validation="internal")
summary(a)#2 optimises
clusters_H <- cutree(clusters_H,2)


library(SLICER)
library(lle)
traj = t(simData)
genes <- select_genes(traj)
traj <- traj[,genes]
k = select_k(traj, kmin=5,kmax=50)
traj_lle = lle(traj, m=2, k)$Y
traj_graph = conn_knn_graph(traj_lle,5)
ends = find_extreme_cells(traj_graph, traj_lle)
start = 7
cells_orderedSLICER = cell_order(traj_graph, start)
pos <- 1:60
pseudoTimesSLICER <- sort(process_distance(traj_graph, start))
assign_branches(traj_graph,start=7)
simSLICER <- simData[,cells_orderedSLICER]
#SLICER thinks there are several branches
write.table(simSLICER,file = "simSLICER2.csv",sep=",",row.names=F,col.names=F)
write.table(pseudoTimesSLICER,file = "pseudoTimesSimSLICER2.csv",sep=",",row.names=F,col.names=F)
###clusters from ordering with SLICER and GPClust
GP_Clust <- unlist(read.table("GPClustSim2Slicer.csv",sep=",",header=F))




#ordering with DeLorean
#uncomment, option (below) to load file with cluster allocations
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
#   sigma.tau=0.5,
#   length.scale=1.5,
#   model.name='exact-sizes')
# simDL <- fit.dl(dl, method='vb')
# simDL<- examine.convergence(simDL)
# plot(simDL, type='Rhat')
# sample.iter=simDL$best.sample
# simDL$samples.l$tau
# pseudoTimesDL <- simDL$samples.l$tau$tau[simDL$samples.l$tau$iter == sample.iter]
# pseudoTimeOrder <- sort(pseudoTimesDL,index.return=T)
# write.table(pseudoTimeOrder$x,file="deLoreanTSim2.csv", sep=",",row.names=F,col.names=F)
# write.table(simData[,pseudoTimeOrder$ix],file="deLoreanSim2.csv", sep=",",row.names=F,col.names=F)
#

GP_Clust1 <- unlist(read.table("GPClustSim2DL.csv",sep=",",header=F))


# ordering and clustering with monocle2

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
cds <- newCellDataSet(as.matrix(simData), phenoData = pd, featureData = fd,expressionFamily = gaussianff())
cds <-  reduceDimension(cds,pseudo_expr = 0,verbose = T,norm_method = "none")
cds <- orderCells(cds)
full_model_fits <- fitModel(cds)
expression_curve_matrix <- responseMatrix(full_model_fits)

distMonocle <- as.dist(1-cor(t(expression_curve_matrix))/2)
sil <- c()
clusters_Monocle<-list()
for (j in 2:15){
  clusters_Monocle[[j-1]] <-  clusterGenes(expression_curve_matrix,k=j)$cluster
  sil = c(sil,mean(silhouette(clusters_Monocle[[j-1]],distMonocle)[,3]))}
clusters_Monocle = clusters_Monocle[[which.max(sil)]]


#import clusters found by SIMLR
clustersSIMLR <- unlist(read.table("sim2_SIMLR.csv"))
#load GPseudoClust summary clusterings
GPseudoClust_lmkk <- read.table('sim2lmkk.csv')
load("npSumClustSim2.rda")

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
write.table(evalMatSim,file="clusterCompSim2.csv",sep=",",row.names=F)
