function [elapsedTime] = compTimeSimLarge(fileNameA,identifiers)
%This function runs the GPSeudoClust algorithm 1 time each on each of
%the 10 simulated larger data sets. It runs only 1 subsampled chain per data set
%and number of subsampled cells (this is therefore for timing only, to analyse
%data run a recommended minimum of 12 chains in parallel. 

%We run this function on different numbers of cores: 1,2,3,4,5,6, as
%multi-threading of matrix computations allows speed-up of runtimes. 

addpath(genpath('../'));
fHandle          =  @GPseudoClust2;  
nSamples         = 5000;  
verbose          = false;    
permuteData      = true;
inputSeed   = NaN;
b           = 0.01;
elapsedTime = zeros(10,5);
adjustForCellSize = false;
for j = 1:10
    fileName         = sprintf('simLarge%d.csv',j);
    A = csvread(fileName);
    [nGenes nCells] = size(A);
    a = nCells/3;
    for k=1:5
        captureTimes     = [ones([1,10*k]),repmat(2,[1,10*k]),repmat(3,[1,10*k])];
        ab = [randsample(1:a,10*k),randsample((a+1):(2*a),10*k),randsample((2*a+1):(3*a),10*k)];
        tic
        feval(fHandle, fileName, identifiers(j), nSamples, verbose,inputSeed,...
                    permuteData,captureTimes,b,adjustForCellSize,ab);
        elapsedTime(j,k) = toc;
    end
end
save(fileNameA,'elapsedTime');
end
