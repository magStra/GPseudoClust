%rerun for timing
addpath(genpath('../'));
captureTimes     = [repmat(2,[1,15]),repmat(4,[1 15]),repmat(6,[1 15])];
fHandle          = @GPseudoClust2;  
fileName         = 'Shalek.csv';
uniqueIdentifiers = 100:112;
nSamples         = 4000; %same as for untimed
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
runTimes = zeros(1,12);
parfor j = 1:12
    %each chain on a single core
subS = [randsample(125:189,15),randsample(190:249,15),randsample(250:307,15)];
%cells from capture times 2h, 4h, and 6h; we run each chain on a different random subset of the cells
tic
feval(fHandle, fileName, uniqueIdentifiers(j), nSamples,verbose, inputSeed,permuteData,captureTimes,...
b,adjustForCellSize,subS);
runTimes(j) = toc;
end
save('runTimesShalek.mat','runTimes');
median(runTimes)/60