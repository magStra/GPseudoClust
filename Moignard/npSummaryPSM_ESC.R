#Computing summary PSMs

PSMs_ESC1 <- readMat("PSMs_ESC_1.mat")
PSMs_ESC1 <- PSMs_ESC1$PSMs.ESC.1
npSummaryPSM1 <- processPSMs(PSMs_ESC1)
writeMat('npSummaryPSM_ESC1.mat',npSummaryPSM1=npSummaryPSM1)

PSMs_ESC2 <- readMat("PSMs_ESC_2.mat")
PSMs_ESC2 <- PSMs_ESC2$PSMs.ESC.1
npSummaryPSM2 <- processPSMs(PSMs_ESC2)
writeMat('npSummaryPSM_ESC2.mat',npSummaryPSM2=npSummaryPSM2)

PSMs_ESC3 <- readMat("PSMs_ESC_3.mat")
PSMs_ESC3 <- PSMs_ESC3$PSMs.ESC.3
npSummaryPSM3 <- processPSMs(PSMs_ESC3)
writeMat('npSummaryPSM_ESC3.mat',npSummaryPSM3=npSummaryPSM3)

#Checking convergence
convESC1 <- checkConvergence(PSMs_ESC1)
save(convESC1,file="convESC1.rda")
convESC2 <- checkConvergence(PSMs_ESC2)
save(convESC2,file="convESC2.rda")
convESC3 <- checkConvergence(PSMs_ESC3)
save(convESC3,file="convESC3.rda")

alphasESC1 <- as.matrix(read.table("alphasESC1.csv",sep=","))
convESC1Alpha4 <- RhatConcentration(alphasESC1,4)
convESC1Alpha8 <- RhatConcentration(alphasESC1,8)
convESC1Alpha12 <- RhatConcentration(alphasESC1,12)
convESC1Alpha24 <- RhatConcentration(alphasESC1,24)
convESC1Alpha32 <- RhatConcentration(alphasESC1,32)
convESC1Alpha48 <- RhatConcentration(alphasESC1,48)

alphasESC2 <- as.matrix(read.table("alphasESC2.csv",sep=","))
convESC2Alpha4 <- RhatConcentration(alphasESC2,4)
convESC2Alpha8 <- RhatConcentration(alphasESC2,8)
convESC2Alpha12 <- RhatConcentration(alphasESC2,12)
convESC2Alpha24 <- RhatConcentration(alphasESC2,24)
convESC2Alpha32 <- RhatConcentration(alphasESC2,32)
convESC2Alpha48 <- RhatConcentration(alphasESC2,48)

alphasESC3 <- as.matrix(read.table("alphasESC3.csv",sep=","))
convESC3Alpha4 <- RhatConcentration(alphasESC3,4)
convESC3Alpha8 <- RhatConcentration(alphasESC3,8)
convESC3Alpha12 <- RhatConcentration(alphasESC3,12)
convESC3Alpha24 <- RhatConcentration(alphasESC3,24)
convESC3Alpha32 <- RhatConcentration(alphasESC3,32)
convESC3Alpha48 <- RhatConcentration(alphasESC3,48)

library(ggplot2)

nChains <- c(rep(4,10),rep(8,10),rep(12,10),rep(24,10),rep(32,10),rep(48,10))
convESC1Alpha <- c(convESC1Alpha4,convESC1Alpha8,convESC1Alpha12,convESC1Alpha24,convESC1Alpha32,convESC1Alpha48)
group <- factor(c(rep("4",10),rep("8",10),rep("12",10),rep("24",10),rep("32",10),rep("48",10)),
                levels = c("4","8","12","24","32","48"))
convESC1AlphaDF <- data.frame(convESC1Alpha=convESC1Alpha,nChains=nChains,group=group)
p <- ggplot(convESC1AlphaDF,aes(x = group,y=convESC1Alpha, group=nChains)) +geom_violin()+geom_jitter(shape=16, position=position_jitter(0.2))+
  xlab("number of chains") + ylab("extended GR statistic")+theme(text=element_text(size=12),axis.text=element_text(size=12))
p

nChains <- c(rep(4,10),rep(8,10),rep(12,10),rep(24,10),rep(32,10),rep(48,10))
convESC2Alpha <- c(convESC2Alpha4,convESC2Alpha8,convESC2Alpha12,convESC2Alpha24,convESC2Alpha32,convESC2Alpha48)
group <- factor(c(rep("4",10),rep("8",10),rep("12",10),rep("24",10),rep("32",10),rep("48",10)),
                levels = c("4","8","12","24","32","48"))
convESC2AlphaDF <- data.frame(convESC2Alpha=convESC2Alpha,nChains=nChains,group=group)
p <- ggplot(convESC2AlphaDF,aes(x = group,y=convESC2Alpha, group=nChains)) +geom_violin()+geom_jitter(shape=16, position=position_jitter(0.2))+
  xlab("number of chains") + ylab("extended GR statistic")+theme(text=element_text(size=12),axis.text=element_text(size=12))
p


nChains <- c(rep(4,10),rep(8,10),rep(12,10),rep(24,10),rep(32,10),rep(48,10))
convESC3Alpha <- c(convESC3Alpha4,convESC3Alpha8,convESC3Alpha12,convESC3Alpha24,convESC3Alpha32,convESC3Alpha48)
group <- factor(c(rep("4",10),rep("8",10),rep("12",10),rep("24",10),rep("32",10),rep("48",10)),
                levels = c("4","8","12","24","32","48"))
convESC3AlphaDF <- data.frame(convESC3Alpha=convESC3Alpha,nChains=nChains,group=group)
p <- ggplot(convESC3AlphaDF,aes(x = group,y=convESC3Alpha, group=nChains)) +geom_violin()+geom_jitter(shape=16, position=position_jitter(0.2))+
  xlab("number of chains") + ylab("extended GR statistic")+theme(text=element_text(size=12),axis.text=element_text(size=12))
p

load("convESC1.rda")
dist_Frobenius <- c(convESC1$distPY,convESC1$distDP)/nrow(convESC1$results[[1]]$weightedPSM)^2
method <- c(rep("DPM+PEAR",length(convESC1$distPY)),rep("PY+PEAR",length(convESC1$distPY)))
number_Chains <- seq(3,96,1)
dist_Frobenius <- data.frame(dist_Frobenius=dist_Frobenius,number_Chains=number_Chains,method=method)
p <- ggplot(dist_Frobenius,aes(x=number_Chains,y=dist_Frobenius)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "Frobenius distance/# matrix elements",breaks = c(0.0001,0.0002,0.0003,0.0004,0.0005,0.001,0.0015))+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(3, 96),breaks=c(3,seq(12,96,12)))+theme(legend.position="top")
p

cophPY <- convESC1$coph[[length(convESC1$coph)]]
cophDPM <- convESC1$coph_DP[[length(convESC1$coph)]]
method <- c(rep("DPM+PEAR",length(cophPY)),rep("PY+PEAR",length(cophPY)))
coph <- c(cophPY,cophDPM)
number_Chains <- seq(2,96,1)
cophDF <- data.frame(coph=coph,method=method,number_Chains=number_Chains)
p <- ggplot(cophDF,aes(x=number_Chains,y=coph)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "cophenetic correlation")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(2, 96),breaks=c(2,seq(12,96,12)))+theme(legend.position="top")
p

cophPYRel <- cophPY/cophPY[length(cophPY)]
cophDPMRel <- cophDPM/cophDPM[length(cophDPM)]
coph <- c(cophPYRel,cophDPMRel)
number_Chains <- seq(2,96,1)
cophDF <- data.frame(coph=coph,method=method,number_Chains=number_Chains)
p <- ggplot(cophDF,aes(x=number_Chains,y=coph)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "ratio of cophenetic correlation")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(2, 96),breaks=c(2,seq(12,96,12)))+theme(legend.position="top")
p


load("convESC2.rda")
dist_Frobenius <- c(convESC2$distPY,convESC2$distDP)/nrow(convESC2$results[[1]]$weightedPSM)^2
method <- c(rep("DPM+PEAR",length(convESC2$distPY)),rep("PY+PEAR",length(convESC2$distPY)))
number_Chains <- seq(3,96,1)
dist_Frobenius <- data.frame(dist_Frobenius=dist_Frobenius,number_Chains=number_Chains,method=method)
p <- ggplot(dist_Frobenius,aes(x=number_Chains,y=dist_Frobenius)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "Frobenius distance/# matrix elements")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(3, 96),breaks=c(3,seq(12,96,12)))+theme(legend.position="top")
p

cophPY <- convESC2$coph[[length(convESC2$coph)]]
cophDPM <- convESC2$coph_DP[[length(convESC2$coph)]]
method <- c(rep("DPM+PEAR",length(cophPY)),rep("PY+PEAR",length(cophPY)))
coph <- c(cophPY,cophDPM)
number_Chains <- seq(2,96,1)
cophDF <- data.frame(coph=coph,method=method,number_Chains=number_Chains)
p <- ggplot(cophDF,aes(x=number_Chains,y=coph)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "cophenetic correlation")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(2, 96),breaks=c(2,seq(12,96,12)))+theme(legend.position="top")
p

cophPYRel <- cophPY/cophPY[length(cophPY)]
cophDPMRel <- cophDPM/cophDPM[length(cophDPM)]
coph <- c(cophPYRel,cophDPMRel)
number_Chains <- seq(2,96,1)
cophDF <- data.frame(coph=coph,method=method,number_Chains=number_Chains)
p <- ggplot(cophDF,aes(x=number_Chains,y=coph)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "ratio of cophenetic correlation")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(2, 96),breaks=c(2,seq(12,96,12)))+theme(legend.position="top")
p


load("convESC3.rda")
dist_Frobenius <- c(convESC3$distPY,convESC3$distDP)/nrow(convESC3$results[[1]]$weightedPSM)^2
method <- c(rep("DPM+PEAR",length(convESC3$distPY)),rep("PY+PEAR",length(convESC3$distPY)))
number_Chains <- seq(3,96,1)
dist_Frobenius <- data.frame(dist_Frobenius=dist_Frobenius,number_Chains=number_Chains,method=method)
p <- ggplot(dist_Frobenius,aes(x=number_Chains,y=dist_Frobenius)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "Frobenius distance/# matrix elements",breaks = c(0.0001,0.0002,0.0003,0.0004,0.0005,0.001,0.0015))+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(3, 96),breaks=c(3,seq(12,96,12)))+theme(legend.position="top")
p

cophPY <- convESC3$coph[[length(convESC3$coph)]]
cophDPM <- convESC3$coph_DP[[length(convESC3$coph)]]
method <- c(rep("DPM+PEAR",length(cophPY)),rep("PY+PEAR",length(cophPY)))
coph <- c(cophPY,cophDPM)
number_Chains <- seq(2,96,1)
cophDF <- data.frame(coph=coph,method=method,number_Chains=number_Chains)
p <- ggplot(cophDF,aes(x=number_Chains,y=coph)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "cophenetic correlation")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(2, 96),breaks=c(2,seq(12,96,12)))+theme(legend.position="top")
p

cophPYRel <- cophPY/cophPY[length(cophPY)]
cophDPMRel <- cophDPM/cophDPM[length(cophDPM)]
coph <- c(cophPYRel,cophDPMRel)
number_Chains <- seq(2,96,1)
cophDF <- data.frame(coph=coph,method=method,number_Chains=number_Chains)
p <- ggplot(cophDF,aes(x=number_Chains,y=coph)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "ratio of cophenetic correlation")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(2, 96),breaks=c(2,seq(12,96,12)))+theme(legend.position="top")
p
