addpath(genpath('../'));
parpool(12);
set(0,'defaultFigureUnits','centimeters');
parfor k = 1:24
%for k=1:1
[mu,nClusters,b,c,sigma2,z,nGenes] = generateDataLinear(sprintf('simLinear%d',k));
captureTimes     = [ones([1,10]),repmat(2,[1,10]),repmat(3,[1,10])];
%capture times for the subsampled cells
fHandle          =  @GPseudoClust2;  
fileName         = sprintf('simLinear%d.csv',k);
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
PSM = psm(dlmread(sprintf('simLinear%d_Results_Chain%d.csv',k,j),',',[49 1 99 nGenes]));
csvwrite(sprintf('simLinear%d_Results_PSM_Chain%d.csv',k,j),PSM);
end
toc
%on BSU-CPU, 24 chains sequentially on single core
%Elapsed time is 1029.697931 seconds.
% Elapsed time is 1085.360844 seconds.
% Elapsed time is 1164.067461 seconds.
% Elapsed time is 1220.618844 seconds.
% Elapsed time is 1225.002997 seconds.
% Elapsed time is 1274.030854 seconds.
% Elapsed time is 1345.800029 seconds.
% Elapsed time is 1374.053091 seconds.
% Elapsed time is 1386.383470 seconds.
% Elapsed time is 1415.837329 seconds.
% Elapsed time is 1626.237219 seconds.
% Elapsed time is 1114.406211 seconds.
% Elapsed time is 1144.015567 seconds.
% Elapsed time is 1155.673016 seconds.
% Elapsed time is 1219.841579 seconds.
% Elapsed time is 1043.169583 seconds.
% Elapsed time is 1243.154452 seconds.
% Elapsed time is 1277.231984 seconds.
% Elapsed time is 1337.380909 seconds.
% Elapsed time is 1462.271476 seconds.
% Elapsed time is 1445.149209 seconds.
% Elapsed time is 1588.862357 seconds.
% Elapsed time is 1592.745456 seconds.

end   

%computation times
mean([1029.697931, 1085.360844,1164.067461,1220.618844,1225.002997,1274.030854,... 
    1345.800029,1374.053091,1386.383470,1415.837329 ,1626.237219,1114.406211,1144.015567,... 
1155.673016,1219.841579,1043.169583,1243.154452,1277.231984,1337.380909 ,1462.271476,1445.149209,... 
1588.862357,1592.745456])/60

median([1029.697931, 1085.360844,1164.067461,1220.618844,1225.002997,1274.030854,... 
    1345.800029,1374.053091,1386.383470,1415.837329 ,1626.237219,1114.406211,1144.015567,... 
1155.673016,1219.841579,1043.169583,1243.154452,1277.231984,1337.380909 ,1462.271476,1445.149209,... 
1588.862357,1592.745456])/60

