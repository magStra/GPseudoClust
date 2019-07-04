captureTimes     = [ones([1,10]),repmat(2,[1,10]),repmat(3,[1,10])];
%capture times for the subsampled cells
fHandle          =  @GPseudoClust2;  
fileName         = 'simDataClust.csv';
uniqueIdentifiers    = 1:24;
nSamples             = 1000;  
verbose          = false; 
inputSeed        = NaN;    
permuteData      = true;
b = 0.01;
adjustForCellSize = false;

parpool(12);%adjust to number of cores
parfor j = 1:24
tic
subS = [randsample(1:20,10),randsample(21:40,10),randsample(41:60,10)];
feval(fHandle, fileName, uniqueIdentifiers(j), nSamples, verbose,inputSeed,...
                permuteData,captureTimes,b,adjustForCellSize,subS);
toc
end
A = [];
for j=1:24
    A = [A;dlmread(sprintf('simDataClust_Results_Chain%d.csv',j),',',[99 1 199 52])];
end
PSM_total = psm(A);
csvwrite('PSM_simClustSubset.csv',PSM_total);

% Elapsed time is 136.989183 seconds.
% Elapsed time is 138.100493 seconds.
% Elapsed time is 141.797974 seconds.
% Elapsed time is 141.904483 seconds.
% Elapsed time is 143.452588 seconds.
% Elapsed time is 147.053516 seconds.
% Elapsed time is 132.770193 seconds.
% Elapsed time is 128.764634 seconds.
% Elapsed time is 126.687226 seconds.
% Elapsed time is 132.985354 seconds.
% Elapsed time is 136.275269 seconds.
% Elapsed time is 133.806262 seconds.
% Elapsed time is 132.182365 seconds.
% Elapsed time is 135.210237 seconds.

%compute individual PSMs to use the different methods
PSM = zeros(52,52,24);
for j=1:24
    PSM(:,:,j) = psm(dlmread(sprintf('simDataClust_Results_Chain%d.csv',j),',',[99 1 199 52]));
end
save('PSMSimNoDropout.mat','PSM');
load('PSMSimNoDropout.mat');

clusterSolution = computeSummaryPSM_lmkk(PSM,2:10);
csvwrite('simlmkk.csv',clusterSolution.clustering);
