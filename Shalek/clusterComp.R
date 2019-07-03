#comparing clustering of Shalek data

setwd("~/Documents/pseudoClust/go_analysis")
source("~/Documents/GPseudoClustCode/GPseudoClust/GO_analysis/go_functions.R")
## -- read clusters
#pseudoClust

load("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/GPseudoClustShalek.rda")
load("~/Documents/pseudoClust/go_analysis/clustersKMeans.rda")
load("~/Documents/pseudoClust/go_analysis/clustersH.rda")
load("~/Documents/pseudoClust/go_analysis/GPClusters.rda")
load("~/Documents/pseudoClust/go_analysis/GPClusters1.rda")
load("~/Documents/pseudoClust/go_analysis/clustersMonocle.rda")
load("~/Documents/pseudoClust/go_analysis/OrigShalekClust.rda")
## -- read GO groups
ensemble.go <- read.csv("~/Documents/GPseudoClustCode/GPseudoClust/GO_analysis/ensembl.mouse_aj_v1.go.txt",
                        sep="\t",stringsAsFactors=FALSE)

head(ensemble.go)
colnames(ensemble.go) <- c("GO","gene")
ensemble.go$"gene" <- toupper(ensemble.go$"gene")

load("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/GPseudoClustShalek.rda")
PseudoClust.ccg <- full.cluster.eval(clustersGPseudo,ensemble.go,F)
PseudoClust.ccg$prop.significant.adj

load("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/ClustShalek_SIMLR.rda")
SIMLRClust.ccg <- full.cluster.eval(clustersSIMLR,ensemble.go,F)
SIMLRClust.ccg$prop.significant.adj

Mclust.ccg <- full.cluster.eval(clustersMclust,ensemble.go,F)
Mclust.ccg$prop.significant.adj

KMeans.ccg <- full.cluster.eval(clustersKMeans,ensemble.go,F)
KMeans.ccg$prop.significant.adj

H.ccg <- full.cluster.eval(clustersH,ensemble.go,F)
H.ccg$prop.significant.adj

GPdeL.ccg <- full.cluster.eval(GPclusters,ensemble.go,F)
GPdeL.ccg$prop.significant.adj

GPSL.ccg <- full.cluster.eval(GPclusters1,ensemble.go,F)
GPSL.ccg$prop.significant.adj

Mon.ccg <- full.cluster.eval(clustersMonocle,ensemble.go,F)
Mon.ccg$prop.significant.adj


O.ccg <- full.cluster.eval(OrigShalekClust,ensemble.go,F)
O.ccg$prop.significant.adj

#boostrapping
Pseudo <- c()
SIMLR <- c()
MCl <- c()
PAM <- c()
H <-c()
GPdeL <- c()
GPSL <- c()
Mon <-c()
O <- c()
for (j in 1:100)
{
  print(j)
  PseudoClust.ccg <- full.cluster.eval(clustersGPseudo,ensemble.go,F)
  Pseudo <- c(Pseudo,PseudoClust.ccg$prop.significant.adj)
  SIMLRClust.ccg <- full.cluster.eval(clustersSIMLR,ensemble.go,F)
  SIMLR <- c(SIMLR,SIMLRClust.ccg$prop.significant.adj)
  Mclust.ccg <- full.cluster.eval(clustersMclust,ensemble.go,F)
  MCl <- c(MCl,Mclust.ccg$prop.significant.adj)
  KMeans.ccg <- full.cluster.eval(clustersKMeans,ensemble.go,F)
  PAM <- c(PAM,KMeans.ccg$prop.significant.adj)
  H.ccg <- full.cluster.eval(clustersH,ensemble.go,F)
  H <- c(H,H.ccg$prop.significant.adj)
  GPdeL.ccg <- full.cluster.eval(GPclusters,ensemble.go,F)
  GPdeL <- c(GPdeL,GPdeL.ccg$prop.significant.adj)
  GPSL.ccg <- full.cluster.eval(GPclusters1,ensemble.go,F)
  GPSL <- c(GPSL,GPSL.ccg$prop.significant.adj)
  Mon.ccg <- full.cluster.eval(clustersMonocle,ensemble.go,F)
  Mon <- c(Mon,Mon.ccg$prop.significant.adj)
  O.ccg <- full.cluster.eval(OrigShalekClust,ensemble.go,F)
  O <- c(O,O.ccg$prop.significant.adj)
}


saveMat <- cbind(O,MCl,SIMLR,PAM,H,GPSL,GPdeL,Mon,Pseudo)
write.table(saveMat,"~/Documents/GPseudoClustCode/comparisons/Shalek_GO.csv",sep=",",row.names=F,col.names=F)

#compute ARI between the different clusterings
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

listOfAllocations = list(clusters_mclust,allocsSIMLR,clusters_KMeans,clusters_H,GP_Clust,GP_Clust1,clusters_Monocle,
                         allocsPseudo)
saveMat1 <- pairwiseARI(listOfAllocations)
write.table(saveMat1,"~/Documents/GPseudoClustCode/comparisons/Shalek_ARI.csv",sep=",",row.names=F,col.names=F)
