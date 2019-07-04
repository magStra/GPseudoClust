%E14tg2a cell line
addpath(genpath(pwd));
captureTimes = csvread('StumpfCT_E.csv');
fHandle          = @GPseudoClust2;  
fileName         = 'Stumpf_E.csv';
nSamples         = 4000; % Number of thinned samples
verbose          = false; % Whether or not to print output to screen
inputSeed        = NaN;     % Set to NaN if you do not want a seed, this will allocate a seed
permuteData      = true;%if true, the sampler starts from a randomly permuted order
b                = 0.01;
adjustForCellSize = false;
CT = captureTimes;
uCT     = unique(CT);
CT_a    = [];
A = {};
for k = 1:length(uCT)
    CT_a = [CT_a;repmat(uCT(k),8,1)];
end
xx = 1:length(CT);
PSMs = zeros(94,94,96);
parpool(12);
tic
parfor j = 1:96
    subS = [];
    for k = 1:length(uCT)
        subS = [subS,randsample(xx(CT==uCT(k)),8)];
    end
    feval(fHandle, fileName, j, nSamples, verbose,inputSeed,permuteData,CT_a,...
                b,adjustForCellSize,subS);
     A = dlmread(sprintf('Stumpf_E_Results_Chain%d.csv',j),',',[1999 1 3999 94]); 
    PSMs(:,:,j) = psm(A);       
end
toc
save('PSMs_StumpfE','PSMs');

