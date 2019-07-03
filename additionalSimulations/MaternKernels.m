addpath(genpath('../'));
set(0,'defaultFigureUnits','centimeters');
%parpool(12);
parfor k = 1:24
%for k=1:1
[mu,nClusters,sigmaW2,sigmaE2,L,l,sigmaw2,nGenes] = generateDataMatern(sprintf('simMatern%d',k));
captureTimes     = [ones([1,10]),repmat(2,[1,10]),repmat(3,[1,10])];
%capture times for the subsampled cells
fHandle          =  @GPseudoClust2;  
fileName         = sprintf('simMatern%d.csv',k);
uniqueIdentifiers    = 1:24;
nSamples             = 1000;  
verbose          = false; 
inputSeed        = NaN;    
permuteData      = true;
b = 0.01;
adjustForCellSize = false;

tic
for j = 1:24
subS = [randsample(1:20,10),randsample(21:40,10),randsample(41:60,10)];
feval(fHandle, fileName, uniqueIdentifiers(j), nSamples, verbose,inputSeed,...
                permuteData,captureTimes,b,adjustForCellSize,subS);
PSM = psm(dlmread(sprintf('simMatern%d_Results_Chain%d.csv',k,j),',',[49 1 99 nGenes]));
csvwrite(sprintf('simMatern%d_Results_PSM_Chain%d.csv',k,j),PSM);
end
toc
end   
%on BSU-CPU - one core per data set for the 24 chains altogether.
% Elapsed time is 971.406876 seconds.
% Elapsed time is 1063.976604 seconds.
% Elapsed time is 1085.980943 seconds.
% Elapsed time is 1159.504072 seconds.
% Elapsed time is 1229.129204 seconds.
% Elapsed time is 1380.327837 seconds.
% Elapsed time is 1384.666488 seconds.
% Elapsed time is 1415.365555 seconds.
% Elapsed time is 1420.310998 seconds.
% Elapsed time is 1504.627993 seconds.
% Elapsed time is 1548.163702 seconds.
% Elapsed time is 1607.060753 seconds.
% Elapsed time is 1246.039825 seconds.
% Elapsed time is 1362.085246 seconds.
% Elapsed time is 1193.159473 seconds.
% Elapsed time is 1149.053228 seconds.
% Elapsed time is 1322.513325 seconds.
% Elapsed time is 1125.651437 seconds.
% Elapsed time is 1245.050740 seconds.
% Elapsed time is 1309.836457 seconds.
% Elapsed time is 1310.230701 seconds.
% Elapsed time is 1280.270415 seconds.
% Elapsed time is 1289.873336 seconds.
% Elapsed time is 1575.387019 seconds.

computeTimes = [971.406876,1063.976604,1085.980943,1159.504072,1229.129204,...
1380.327837,1384.666488 ,1415.365555,1420.310998,1504.627993 ,1548.163702 ...
,1607.060753 ,1246.039825 ,1362.085246,1193.159473 ,1149.053228,1322.513325 ...
,1125.651437 ,1245.050740 ,1309.836457,1310.230701,1280.270415,1289.873336,1575.387019];

mean(computeTimes)
median(computeTimes)
