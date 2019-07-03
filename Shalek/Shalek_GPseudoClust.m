addpath(genpath('../'));
captureTimes     = [repmat(2,[1,15]),repmat(4,[1 15]),repmat(6,[1 15])];
fHandle          = @GPseudoClust2;  
fileName         = 'Shalek.csv';
uniqueIdentifiers = 1:96;
nSamples         = 4000; %number of thinned samples
verbose          = false; % Whether or not to print output to screen
inputSeed        = NaN;     % Set to NaN if you do not want a seed, this will allocate a seed
%determined by the clock and a unique deviation for each chain depending on
%the uniqueIdentifier of that chain
permuteData      = true;%if true, the sampler starts from a randomly permuted order
%(whether this is set to true or false is irrelevant here, as we use randomly
%chosen unordered subsamples of cells)
b = 0.01;%recommended standard setting
adjustForCellSize = false;% set to true if the data set has not been adjusted for 
%cell size and you want it to be adjusted (our data have already been
%adjusted)
parpool(12);%adjust to the number of cores
%If you use 12 cores and set parpool(12), then each chain will run on one
%core. For k*12 cores uses maxNumCompThreads(k) within the parfor loop to
%for speed increase using multithreading for matrix computations. 
tic
parfor j = 1:96
    %maxNumCompThreads(k);
subS = [randsample(125:189,15),randsample(190:249,15),randsample(250:307,15)];
%cells from capture times 2h, 4h, and 6h; we run each chain on a different random subset of the cells
feval(fHandle, fileName, uniqueIdentifiers(j), nSamples,verbose, inputSeed,permuteData,captureTimes,...
b,adjustForCellSize,subS);
end
toc
%PSMs for each chain separately
PSMsShalek = zeros(74,74,96);
for j = 1:96
    B = dlmread(sprintf('Shalek_Results_Chain%d.csv',j),',',[2000 1 3999 74]);
    PSMsShalek(:,:,j) = psm(B);
end
save('PSMsShalek.mat','PSMsShalek');

%read in concentration parameter and save as csv file
alphasShalek = zeros(4000,96);
for j = 1:96
    alphasShalek(:,j) = dlmread(sprintf('Shalek_Results_Chain%d.csv',j),',',[0 0 3999 0]);
end
csvwrite('alphasShalek.csv',alphasShalek);

%for the postprocessing of the PSMs first run R code 
%ShalekNonparametricSummaryPSM.R, which uses our R package
%nonparametricSummaryPSM

%then run Shalek_postprocess_plot

