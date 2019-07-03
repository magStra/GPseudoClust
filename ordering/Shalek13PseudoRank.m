%Sampling orders for the Shalek13 data set using GPseudoRank
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 14)
set(0,'defaultfigurecolor',[1 1 1]);
addpath(genpath(pwd));
fHandle          = @pseudoRank; 
fileName         = 'Shalek13PreOrderedAV.csv';
nSamples         = 20000;  
priorLogSDLength = 0.01;%standard deviation of prior of log-length scale, recommended value
paramSamplingFreq =10; %sampling parameters for every kth order
verbose          = true; % Whether or not to print output to screen
initialise       = true;  % If you have to stop the sampler, and then wish to rerun at a later date, set this to false 
thinningFreq     = 1;     % If set to X, records every X-th sample
inputSeed        = NaN;  %use clock seed with chain-dependent offset
delta            = [1/8000,1/8000];
pp               = [0 0 0.4995 0.4995 0.002];
permuteData      = true;
captureTimes     = [repmat(1,1,3),repmat(2,1,3),repmat(3,1,3),...
    repmat(4,1,3),repmat(5,1,3),repmat(6,1,3)];
regInt           = false;%fixed recommended setting
permutationFileName = NaN;
stepSize = [0.1,0.1];%fixed recommended setting
n0 = max(floor(18/4),1);
        n3 = 3;
        n3a = 3;
        kk = 0.01;
        jj = 1;
parpool(3);
for j=1:3
 tic
 feval(fHandle, fileName, j, nSamples, priorLogSDLength, verbose, initialise, thinningFreq, paramSamplingFreq,...
     stepSize,inputSeed, delta, pp,permuteData,captureTimes,regInt,permutationFileName,20000,n0,n3,n3a,jj,kk);
 
 toc
end
 %161 seconds on 1 core if we save at every iteration
 %33 seconds on 1 core if saving at every 10th iteration
 
 
 
 
