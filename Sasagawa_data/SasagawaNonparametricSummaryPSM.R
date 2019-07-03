library(nonparametricSummaryPSM)
library(R.matlab)
library(mcclust)

PSMs_Sasagawa <- readMat("PSMs_mes_sub.mat")
PSMs_Sasagawa <- PSMs_Sasagawa$PSMS
SasagawaPSMs <- processPSMs(PSMs_Sasagawa)
npSummaryPSMSasagawa <- processPSMs(PSMs_Sasagawa)
writeMat('npSummaryPSMSasagawa.mat',npSummaryPSMSasagawa=npSummaryPSMSasagawa)

PSMs_Sasagawa <- readMat("PSMs_mes_sub.mat")
convSasagawa <- checkConvergence(PSMs_Sasagawa$PSMS)
save(convSasagawa,file="convSasagawa.rda")
alphasSasagawa <- as.matrix(read.table("alphasSasagawa.csv",sep=","))
convSasagawaAlpha4 <- RhatConcentration(alphasSasagawa,4)
convSasagawaAlpha8 <- RhatConcentration(alphasSasagawa,8)
convSasagawaAlpha12 <- RhatConcentration(alphasSasagawa,12)
convSasagawaAlpha18 <- RhatConcentration(alphasSasagawa,18)

library(ggplot2)
load("convSasagawa.rda")
dist_Frobenius <- c(convSasagawa$distPY,convSasagawa$distDP)/nrow(npSummaryPSMSasagawa$weightedPSM)^2
method <- c(rep("DPM+PEAR",length(convSasagawa$distPY)),rep("PY+PEAR",length(convSasagawa$distPY)))
number_Chains <- seq(3,36,1)
dist_FrobeniusDF <- data.frame(dist_Frobenius=dist_Frobenius,number_Chains=number_Chains,method=method)
p <- ggplot(dist_FrobeniusDF,aes(x=number_Chains,y=dist_Frobenius)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "Frobenius distance/# matrix elements",breaks = seq(0.00002,0.00012,0.00002))+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(3, 36),breaks=c(3,seq(6,36,6)))+theme(legend.position="top")
p

cophPY <- convSasagawa$coph[[length(convSasagawa$coph)]]
cophDPM <- convSasagawa$coph_DP[[length(convSasagawa$coph)]]
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
convSasagawaAlpha <- c(convSasagawaAlpha4,convSasagawaAlpha8,convSasagawaAlpha12,convSasagawaAlpha18)
group <- factor(c(rep("4",10),rep("8",10),rep("12",10),rep("18",10)),
                levels = c("4","8","12","18"))
convSasagawaAlphaDF <- data.frame(convSasagawaAlpha=convSasagawaAlpha,nChains=nChains,group=group)
p <- ggplot(convSasagawaAlphaDF,aes(x = group,y=convSasagawaAlpha, group=nChains)) +geom_violin()+geom_jitter(shape=16, position=position_jitter(0.2))+
  xlab("number of chains") + ylab("extended GR statistic")+theme(text=element_text(size=12),axis.text=element_text(size=12))
p
