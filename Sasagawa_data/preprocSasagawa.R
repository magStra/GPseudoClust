dataQuartz <- read.csv("~/Documents/GPseudoClustCode/dataQuartz.csv")#change directory
 ###using only cell-cycle related genes
 library(AnnotationDbi)
 library(org.Mm.eg.db)
 #the following lines (6-13) taken from https://github.com/PMBio/scLVM/blob/master/R/scripts/transform_counts_Tcells.R

xxGO <- as.list(org.Mm.egGO2EG)
 cell_cycleEG <-unlist(xxGO['GO:0007049'])
 #get ENSEMBLE ids
 x <- org.Mm.egENSEMBL
 mapped_genes <- mappedkeys(x)
 xxE <- as.list(x[mapped_genes])
 ens_ids_cc<-unlist(xxE[cell_cycleEG])
 cc_gene_indices <- na.omit(match(ens_ids_cc, dataQuartz[,1]))
 ##
 mesCC <- log2(dataQuartz[cc_gene_indices,2:36]+1)
 CC_geneNames <- dataQuartz[cc_gene_indices,1]
 write.table(mesCC,"~/Documents/GPseudoClustCode/mesCC.csv",sep=",",col.names=F,row.names=F)
 write.table(CC_geneNames,"~/Documents/GPseudoClustCode/mesGenesCC.csv",sep=",")

