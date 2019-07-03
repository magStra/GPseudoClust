# #mclust
ShalekData <- read.csv("~/Documents/pseudoRank/Shalek.csv")
ShalekGeneNames <-(ShalekData[,1])
ShalekData <- ShalekData[,2:308]
ShalekData <- ShalekData[,125:307]#only those we also use with GPseudoClust
library(mclust)
library(cluster)
# 
clusters_mclust <- Mclust(ShalekData,G=4:20)$classification
clustersMclust <- list()
for (j in 1:max(clusters_mclust)){
  clustersMclust[[j]] <- ShalekGeneNames[which(clusters_mclust == j)]}

#k-means
library(fpc)
clusters_KMeans <- pamk(ShalekData,4:15)$pamobject$clustering
clusters_KMeans12 <- pam(ShalekData,12)$cluster
clustersKMeans <- list()
for (j in 1:max(clusters_KMeans)){
  clustersKMeans[[j]] <- ShalekGeneNames[which(clusters_KMeans == j)]}
clustersKMeans12 <- list()
for (j in 1:max(clusters_KMeans12)){
  clustersKMeans12[[j]] <- ShalekGeneNames[which(clusters_KMeans12 == j)]}

clusters_H <- hclust(dist(ShalekData))
plot(clusters_H)
clusters_H12 <- cutree(clusters_H,k=12)#to have same number of clusters as with GPseudoRank with subsampling
clusters_H <- cutree(clusters_H,k=4)#from silhouette width with clValid
clustersH <- list()
for (j in 1:max(clusters_H)){
  clustersH[[j]] <- ShalekGeneNames[which(clusters_H == j)]}
clustersH12 <- list()
for (j in 1:max(clusters_H12)){
  clustersH12[[j]] <- ShalekGeneNames[which(clusters_H12 == j)]}
# ########
# 
# 
# library(SLICER)
# library(lle)
# traj = t(ShalekData)
# genes <- select_genes(traj)
# traj <- traj[,genes]
# k = select_k(traj, kmin=5,kmax=50)
# traj_lle = lle(traj, m=2, k)$Y
# traj_graph = conn_knn_graph(traj_lle,5)
# ends = find_extreme_cells(traj_graph, traj_lle)
# start = ends[2]
# cells_orderedSLICER = cell_order(traj_graph, start)
# pseudoTimesSLICER <- sort(process_distance(traj_graph, start))
# assign_branches(traj_graph,start=ends[2])
#                 ###SLICER recognises there is no branching
# ShalekSLICER <- ShalekData[,cells_orderedSLICER]
# # write.csv(ShalekSLICER,file = "~/Documents/pseudoClust/ShalekSLICER.csv")
# # write.csv(pseudoTimesSLICER,file = "~/Documents/pseudoClust/pseudoTimesSLICER.csv")
# # write.csv(shQuote(ShalekGeneNames),file = "~/Documents/pseudoClust/ShalekGeneNames.csv")
# 
# ###clusters from ordering with SLICER and Hensman method
# #ordering and GPclust
# GP_Clust <- c(3.,  2. , 4. , 2. , 3. , 1. , 1.  ,5. , 4. , 1. , 4.,  1. , 1. , 2.,  5. , 1.,  1.,  1.,
# 2.,  1.,  1.,  4. , 2.,  5. , 2. , 3. , 3.,  3. , 3.,  2. , 1.,  1. , 2. , 3. , 1.,  3.,
# 3, 1,  1,  6,  2, 2,  1,  1,  1,  2,  1,  2,  4,  1,  3,  1,  4,  1,
# 1,  6,  3,  2,  4,  5,  1,  5,  3,  1,  4,  1,  2,  1,  2,  1,  1,  3,
# 3,  1) #from Python
# 
# GPclusters <- list()
# for (j in 1:max(GP_Clust)){
#   GPclusters[[j]] <- ShalekGeneNames[which(GP_Clust == j)]}
# 
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
#write.csv(pseudoTimeOrder$x,file="~/Documents/pseudoClust/deLoreanTShalek.csv")
#write.csv(ShalekData[,pseudoTimeOrder$ix],file="~/Documents/pseudoClust/deLoreanShalek.csv")

GP_Clust1 <- c(2,  1,  1,  1,  6,  1,  1,  2,  1,  1,  1,  1,  1,  1,  4,  1,  1,  1,
               1,  1,  1,  1,  1,  2,  1,  5,  3,  2,  5,  1,  1,  1,  1,  3,  1,  3,
               2,  1,  1,  2,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  6,  1,  1,  1,
               1,  4,  5,  1,  2,  2,  1,  4,  3,  1,  1,  1,  1,  1,  1,  1,  1,  2,
               3,  1)
GPclusters1 <- list()
for (j in 1:max(GP_Clust1)){
  GPclusters1[[j]] <- ShalekGeneNames[which(GP_Clust1 == j)]}

# ordering and clustering with monocle2
library(monocle)
load("~/Documents/PseudoClust/ShalekDeLorean.rda")
xx <- c()
yy <- c()
for (j in 1:74){#if we need the same order of the cells and genes, we can't use %in%
  xx <- c(xx,which(rownames(shalek.A.expr) == ShalekGeneNames[j]))}
for (j in 1:dim(ShalekData)[2]){
  yy <- c(yy,which(colnames(shalek.A.expr) == colnames(ShalekData[j])))}
ShalekSubset <- shalek.A.expr[xx,yy]#this is not adjusted

pd <- shalek.A.cell.meta[yy,]
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
clustersMonocle <- list()
for (j in 1:max(clusters_Monocle)){
  clustersMonocle[[j]] <- ShalekGeneNames[which(clusters_Monocle == j)]}
clustersMonocle12 <- list()
clusters_Monocle12 <- clusterGenes(expression_curve_matrix,k=12)$cluster
for (j in 1:max(clusters_Monocle12)){
  clustersMonocle12[[j]] <- ShalekGeneNames[which(clusters_Monocle12 == j)]}

#Shalek data, finding the original gene modules
options(stringsAsFactors = FALSE)
load("~/Documents/pseudoClust/ShalekDeLorean.rda")
ShalekExtract = read.csv("~/Documents/pseudoRank/Shalek.csv")
extractGenes = ShalekExtract$X
A <- (shalek.A.gene.meta$gene)
extractGeneModules <- c()
for (j in 1:length(extractGenes)){
  extractGeneModules[j] = shalek.A.gene.meta$cluster[A == extractGenes[j]]
}
uu <- unique(extractGeneModules)
OrigShalekClust <- list()
for (j in 1:4)
{
  OrigShalekClust[[j]] <- extractGenes[extractGeneModules==uu[j]]
}
save(OrigShalekClust,file="~/Documents/pseudoClust/go_analysis/OrigShalekClust.rda")
save(clustersMclust,clusters_mclust, file="~/Documents/pseudoClust/go_analysis/clusterMClustShalek74.rda")
save(clustersKMeans,clusters_KMeans, file="~/Documents/pseudoClust/go_analysis/clustersKmeans.rda")
save(clustersKMeans12,clusters_KMeans12, file="~/Documents/pseudoClust/go_analysis/clustersKmeans12.rda")
save(clustersH,clusters_H, file="~/Documents/pseudoClust/go_analysis/clustersH.rda")
save(clustersH12,clusters_H12, file="~/Documents/pseudoClust/go_analysis/clustersH12.rda")
save(GPclusters,GP_Clust,file="~/Documents/pseudoClust/go_analysis/GPClusters.rda")
save(GPclusters1,GP_Clust1,file="~/Documents/pseudoClust/go_analysis/GPClusters1.rda")
save(clustersMonocle,clusters_Monocle,file="~/Documents/pseudoClust/go_analysis/clustersMonocle.rda")
save(clustersMonocle12,clusters_Monocle12,file="~/Documents/pseudoClust/go_analysis/clustersMonocle12.rda")
