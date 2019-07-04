library(nonparametricSummaryPSM)
library(R.matlab)
library(mcclust)

PSMs_Stumpf <- array(dim=c(94,94,192))
PSMs_StumpfE <- readMat("PSMs_StumpfE.mat")
PSMs_StumpfR <- readMat("PSMs_StumpfR.mat")
PSMs_Stumpf[,,1:96] <- PSMs_StumpfE$PSMs
PSMs_Stumpf[,,97:192] <- PSMs_StumpfR$PSMs
npSummaryPSMStumpf <- processPSMs(PSMs_Stumpf)
writeMat('npSummaryPSMStumpf.mat',npSummaryPSMStumpf=npSummaryPSMStumpf)



