function [] = simClustManyGenes(nGenes,nSamples)
    addpath(genpath('../'));
    uniqueIdentifier = nGenes;
    inputSeed = sum(100*clock) + sum(100*double(uniqueIdentifier)); 
    simManyGenes(sprintf('simMoreGenes%d',uniqueIdentifier),inputSeed,nGenes);
    fHandle          =  @GPseudoClust2;  
    fileName         = sprintf('simMoreGenes%d.csv',uniqueIdentifier);
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
    parpool(10);
    elapsedTime = zeros(10,1);
    parfor jj = 1:10
         tic
         ab = [randsample(1:a,10),randsample((a+1):(2*a),10),randsample((2*a+1):(3*a),10)];
         feval(fHandle, fileName, jj, nSamples, verbose,inputSeed,...
                        permuteData,captureTimes,b,adjustForCellSize,ab);
        elapsedTime(jj) = toc
        PSMs(:,:,jj) = psm(dlmread(sprintf('simMoreGenes%d_Results_Chain%d.csv',uniqueIdentifier,jj),',',[floor(nSamples/2) 1 nSamples-1 nGenes]));
    end
    pool = gcp('nocreate');
    delete(pool);
    save(sprintf('PSMsMoreGenes%d_%d.mat',nGenes,uniqueIdentifier),'PSMs','elapsedTime');
end



