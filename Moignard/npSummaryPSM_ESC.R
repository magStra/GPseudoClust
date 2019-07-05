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
