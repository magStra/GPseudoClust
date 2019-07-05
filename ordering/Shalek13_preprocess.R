#Data downloaded from http://imlspenticton.uzh.ch:3838/conquer/
#Charlotte, S and Robinson MD (2018). Bias, robustness and scalability in single-cell
#differential expression analysis. Nature Methods 15, pages 255–261.

#For comparability between GPseudoClust and GPseudoRank, we use the same preprocessing steps as in
#Strauss et al. 2018.
Shalek13 <- readRDS("[path to file]/GSE41265.rds")
all <- assay(Shalek13)

all <- log2(all+1)
inds3 <- apply(all,1,function(x){return(sum(x!=0)>18*0.3)})
all1 <- all[inds3,]
inds <- order(apply(all1,1,var),decreasing=T)[1:1000]
inds1 <- order(apply(all1,1,mean),decreasing=T)[1:1000]
inds2 <- intersect(inds,inds1)
Shalek13 <- all1[inds2,]
write.table(rownames(Shalek13),"Shalek13Genes.csv",sep=",",col.names = F,row.names = F)
write.table(as.matrix(Shalek13),"Shalek13.csv",sep=",",col.names = F,row.names = F)

#antiviral score for Shalek13
library(biomaRt)
#core-antiviral genes as in Shalek et al. (2013), Reid and Wernisch (2016), Strauss et al. (2018)
coreAV <- c("TOR3A", "DAXX", "TAP1", "H2-T10", "EIF2AK2","GM4951", "IIGP1", "PYHIN1", "E030037K03RIK",
            "AI607873", "IFI204", "MNDAL", "IFI203", "IFI205", "AW112010", "MS4A4C", "CD274",
            "IFIT2", "IFIT3", "GM14446", "I830012O16RIK", "IFIT1", "IL15RA", "IFIH1", "ZNFX1",
            "ZBP1", "BC006779", "CAR13", "ADAR", "ZUFSP", "FAM26F", "MOV10", "GBP3", "GBP2",
            "IFI44","DDX58","MITD1", "ISG15","CXCL10", "MPA2L", "OASL2", "OASL1", "OAS2",
            "OAS3", "OAS1G","OAS1A", "SAMD9L", "PARP12", "NT5C3", "HERC6","USP18", "CD69", "ETNK1",
            "SLCO3A1", "ISG20", "TRIM30A", "TRIM30D", "GVIN1", "GM8979", "AK172683", "IFITM3", "IRF7", "DDX60", "BST2",
            "IL15", "NLRC5", "CASP11", "STAT1", "PML" ,"UBA7", "TREX1" ,"STAT2" ,"PTTG1" ,"IRGM1", "GM5431" ,
            "IFI47" ,"GM12250", "IRGM2" ,"IGTP", "XAF1", "SLFN5", "SLFN9", "SLFN8" ,"DHX58" ,"IFI35", "RSAD2",
            "CMPK2", "SP140", "SP100", "SERPINA3G", "PHF11", "GM4902", "D14ERTD668E", "RTP4", "DTX3L", "PARP9","USP25", "MX1", "RNASET2A"
)
ensembl = useMart("ensembl")
mouse.ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
cAV <- getBM(attributes = "ensembl_gene_id", filters = 'mgi_symbol',
             values=coreAV,mart=mouse.ensembl)
cAV <- unlist(cAV$ensembl_gene_id)
rowN <- rownames(all)
for (j in 1:length(rowN))
{
  rowN[j] <- unlist(strsplit(cAV[j], split='.', fixed=TRUE))[1]
}
rownames(all) <- rowN
Shalek13AV <- colMeans(all[rownames(all) %in% cAV,])
write.table(Shalek13AV,"Shalek13AV.csv",sep=",",row.names = F,col.names = F)
Shalek13_av <- all[rownames(all) %in% cAV,]
write.table(Shalek13_av,"Shalek13_av.csv",sep=",",row.names = F,col.names = F)
