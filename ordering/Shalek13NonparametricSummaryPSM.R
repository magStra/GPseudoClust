library(nonparametricSummaryPSM)
library(R.matlab)
library(mcclust)

PSMs_Shalek13 <- readMat("Shalek13_PSMs.mat")
PSMs_Shalek13 <- PSMs_Shalek13$PSMs
Shalek13PSMs <- processPSMs(PSMs_Shalek13)
npSummaryPSMShalek13 <- processPSMs(PSMs_Shalek13)
writeMat('npSummaryPSMShalek13.mat',npSummaryPSMShalek13=npSummaryPSMShalek13)

convShalek13 <- checkConvergence(PSMs_Shalek13)
save(convShalek13,file="convShalek13.rda")
alphasShalek13 <- as.matrix(read.table("alphasShalek13.csv",sep=","))
convShalek13Alpha4 <- RhatConcentration(alphasShalek13,4)
convShalek13Alpha8 <- RhatConcentration(alphasShalek13,8)
convShalek13Alpha12 <- RhatConcentration(alphasShalek13,12)
convShalek13Alpha18 <- RhatConcentration(alphasShalek13,18)

library(ggplot2)
load("convShalek13.rda")
dist_Frobenius <- c(convShalek13$distPY,convShalek13$distDP)/nrow(npSummaryPSMShalek13$weightedPSM)^2
method <- c(rep("DPM+PEAR",length(convShalek13$distPY)),rep("PY+PEAR",length(convShalek13$distPY)))
number_Chains <- seq(3,36,1)
dist_FrobeniusDF <- data.frame(dist_Frobenius=dist_Frobenius,number_Chains=number_Chains,method=method)
p <- ggplot(dist_FrobeniusDF,aes(x=number_Chains,y=dist_Frobenius)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "Frobenius distance/# matrix elements",breaks = seq(0.00002,0.00012,0.00002))+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(3, 36),breaks=c(3,seq(6,36,6)))+theme(legend.position="top")
p

cophPY <- convShalek13$coph[[length(convShalek13$coph)]]
cophDPM <- convShalek13$coph_DP[[length(convShalek13$coph)]]
method <- c(rep("DPM+PEAR",length(cophPY)),rep("PY+PEAR",length(cophPY)))
coph <- c(cophPY,cophDPM)
number_Chains <- seq(2,36,1)
cophDF <- data.frame(coph=coph,method=method,number_Chains=number_Chains)
p <- ggplot(cophDF,aes(x=number_Chains,y=coph)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "cophenetic correlation")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(2, 36),breaks=c(2,seq(6,36,6)))+theme(legend.position="top")
p
cophPYRel <- cophPY/cophPY[length(cophPY)]
cophDPMRel <- cophDPM/cophDPM[length(cophDPM)]
coph <- c(cophPYRel,cophDPMRel)
number_Chains <- seq(2,36,1)
cophDF <- data.frame(coph=coph,method=method,number_Chains=number_Chains)
p <- ggplot(cophDF,aes(x=number_Chains,y=coph)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "ratio of cophenetic correlation",breaks=seq(0.94,1.01,0.01))+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(2, 36),breaks=c(2,seq(6,36,6)))+theme(legend.position="top")
p

nChains <- c(rep(4,10),rep(8,10),rep(12,10),rep(18,10))
convShalek13Alpha <- c(convShalek13Alpha4,convShalek13Alpha8,convShalek13Alpha12,convShalek13Alpha18)
group <- factor(c(rep("4",10),rep("8",10),rep("12",10),rep("18",10)),
                levels = c("4","8","12","18"))
convShalek13AlphaDF <- data.frame(convShalek13Alpha=convShalek13Alpha,nChains=nChains,group=group)
p <- ggplot(convShalek13AlphaDF,aes(x = group,y=convShalek13Alpha, group=nChains)) +geom_violin()+geom_jitter(shape=16, position=position_jitter(0.2))+
  xlab("number of chains") + ylab("extended GR statistic")+theme(text=element_text(size=12),axis.text=element_text(size=12))
p


