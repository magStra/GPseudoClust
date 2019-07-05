captureTimes     = [ones([1,10]),repmat(2,[1,10]),repmat(3,[1,10])];
addpath(genpath('../'));
fHandle          =  @GPseudoClust2;  
fileName         = 'simDataClust2.csv';
uniqueIdentifiers    = 1:24;
nSamples             = 1000;  
verbose          = true; 
inputSeed        = NaN;    
permuteData      = true;
b = 0.01;
adjustForCellSize = false;

%parpool(12);%adjust to number of cores
for j = 1:24
tic
subS = [randsample(1:20,10),randsample(21:40,10),randsample(41:60,10)]
disp(subS)
feval(fHandle, fileName, uniqueIdentifiers(j), nSamples, verbose,inputSeed,...
                permuteData,captureTimes,b,adjustForCellSize,subS);
toc
end
A = [];
for j=1:24
    A = [A;dlmread(sprintf('simDataClust2_Results_Chain%d.csv',j),',',[99 1 199 52])];
end
PSM_total = psm(A);
csvwrite('PSM_simClustSubset2.csv',PSM_total);
%compute individual PSMs to use the different methods
PSM = zeros(52,52,24);
for j=1:24
    PSM(:,:,j) = psm(dlmread(sprintf('simDataClust2_Results_Chain%d.csv',j),',',[99 1 199 52]));
end
save('PSMSim2NoDropout.mat','PSM');
load('PSMSim2NoDropout.mat');

clusterSolution = computeSummaryPSM_lmkk(PSM,2:10);
csvwrite('sim2lmkk.csv',clusterSolution.clustering);
