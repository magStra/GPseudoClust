function [] = simClustManyGenes2(uniqueIdentifier)
    %to simulate a dataset with 5000 genes and infer clusters with GPseudoClust
    %first run simManyGenes('simMoreGenes.csv',[choose input seed],5000)
%then run this function with uniqueIdentifiers 1:12 in parallel
    addpath(genpath('../'));
    nGenes = 5000;
    nSamples = 1000;
    inputSeed = sum(100*clock) + sum(100*double(uniqueIdentifier));
    fHandle          =  @GPseudoClust2;  
    fileName         = sprintf('simMoreGenes%d.csv',111);%the same for all subsampled chains
    verbose          = false; 
    inputSeed        = NaN;    
    permuteData      = true;
    b = 0.01;
    adjustForCellSize = false;
    A = csvread(fileName);
    [nGenes nCells] = size(A);
    a = nCells/3;
    PSMs = zeros(nGenes,nGenes,12);
    captureTimes     = [ones([1,10]),repmat(2,[1,10]),repmat(3,[1,10])];
     tic
     ab = [randsample(1:a,10),randsample((a+1):(2*a),10),randsample((2*a+1):(3*a),10)];
     feval(fHandle, fileName, uniqueIdentifier, nSamples, verbose,inputSeed,...
                    permuteData,captureTimes,b,adjustForCellSize,ab);
    elapsedTime = toc
    PSM=psm(dlmread(sprintf('simMoreGenes%d_Results_Chain%d.csv',111,uniqueIdentifier),',',[floor(nSamples/2) 1 nSamples-1 nGenes]));
    save(sprintf('PSMsMoreGenes%d.mat',uniqueIdentifier),'PSM','elapsedTime');
end
