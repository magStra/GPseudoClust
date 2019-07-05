addpath(genpath(pwd))
CT = csvread('ESC_CT_3.csv');
uCT = unique(CT)
captureTimes     = [];
for j = 1:length(uCT)
   a = min(floor(sum(CT==uCT(j))/4),12)
    captureTimes = [captureTimes;repmat(uCT(j),a,1)];
end

fHandle          = @GPseudoClust2;  
fileName         = 'ESC_3.csv';
uniqueIdentifiers = 1:96;
nSamples         = 2000; %number of thinned samples
verbose          = false; % Whether or not to print output to screen
inputSeed        = NaN;     % Set to NaN if you do not want a seed, this will allocate a seed
%determined by the clock and a unique deviation for each chain depending on
%the uniqueIdentifier of that chain
permuteData      = true;%if true, the sampler starts from a randomly permuted order
b = 0.01;
adjustForCellSize = false;% set to true if the data set has not been adjusted for 
%cell size and you want it to be adjusted (our data have already been
%adjusted)
%parpool(16);%adjust to the number of cores
xx = 1:length(CT);
 tic
 parfor j = 1:96
subS = [];
     for k = 1:length(uCT)
         yy = xx(CT == uCT(k));
         a = min(floor(length(yy)/4),12);
         subS = [subS,randsample(yy,a)];
     end
feval(fHandle, fileName, uniqueIdentifiers(j), nSamples,verbose, inputSeed,permuteData,captureTimes,...
b,adjustForCellSize,subS);
end
toc
%PSMs for each chain separately
PSMs_ESC_3 = zeros(42,42,96);
for j = 1:96
    B = dlmread(sprintf('ESC_3_Results_Chain%d.csv',j),',',[1000 1 1999 42]);
    PSMs_ESC_3(:,:,j) = psm(B);
end
save('PSMs_ESC_3.mat','PSMs_ESC_3');
