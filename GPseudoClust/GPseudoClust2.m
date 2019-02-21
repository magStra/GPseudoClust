function [] = GPseudoClust2(fileName, uniqueIdentifier, nSamples,...
verbose, inputSeed, permuteData,captureTimes,b,adjustForCellSize,subSelect)
%This is the main function required to run the GPseudoClust method. 
%input:
%fileName: name of the .csv file containing the mRNA-expression matrix
%uniqueIdentifier: positive integer to identify the chain
%nSamples: number of samples to be drawn from the posterior distribution
%verbose: if true prints acceptance rates
%inputSeed: set for the random number generator, if set to NaN it is
%clock-seeded with a chain-dependent offset
%permuteData: should generally be set to true to make the chain start from
%a random order of the cells
%captureTimes: the capture times of the cells on which the sampler is run
%b: prior variance of log-length scale
%adjustForCellSize: if true the data are adjusted for cell size using a
%method by Anders, S. and Huber,W. (2010) . Differential expression 
%analysis for sequence count data. Genome Biol, 11(10), R106?R106.
%subSelect: the subset of cells on which we run the chain

%%%%%%%
%output files:
%Three .csv files: FileName refers to the fileName without the .csv
%extension:
    %1) FileName_PermutationuniqueIdentifier.csv: contains the permutation
    %of the original order of the cells that was used as the starting order
    %2) FileName_Results_ChainuniqueIdentifier.csv: contains the samples of
    %the cluster allocations, the concentration parameter, and the GP
    %hyperparameters
    %3) FileName_Results_Lk_ChainuniqueIdentifier.csv: the likelihoods of
    %the samples
    
%part of the code is taken or adapted from (see comments in the text
%indicating the lines): 
%Paul DW Kirk. DPMSysBio. GitHub repository. https://github.com/pauldwkirk/DPMSysBio    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clock-seed the random number generator (with a chain-depenedent offset)
if(isnan(inputSeed))
    inputSeed = sum(100*clock) + sum(100*double(uniqueIdentifier)); 
end
rng(inputSeed);

gammaPrior       = [2 4]; 
hyperParameterSamplingFrequency = 10;%sample GP parameters at every 10th 
%iteration 
clusterAllocationFrequency = 5;
thinningFreq = 5;
nSamples = nSamples*5;
stepSize = [0.01,0.1,0.05,0.3];

saveFileName = [strtok(fileName, '.'),'_Results_Chain', num2str(uniqueIdentifier)];
saveFileNameOrder = [strtok(fileName, '.'),'_Results_Orders_Chain', num2str(uniqueIdentifier)];
saveFileNameLk      = [strtok(fileName, '.'),'_Results_Lk_Chain', num2str(uniqueIdentifier)];


% We need to initialise the variables and structures
fHandle = @TimeCourse_pseudotime2;
% Read in data 
allData = importdata(fileName, ',');
try
    data    = allData.data;
catch
    data    = allData;%if there are no gene names and cell names
end
data = data(:,subSelect);
if adjustForCellSize == true
    uniqueCaptureTimes = unique(captureTimes);
    lu = length(uniqueCaptureTimes);
    Data1 = [];
    for j = 1:lu
        aa = data(:,captureTimes==uniqueCaptureTimes(j));
        aaMod = bsxfun(@minus, aa,mean(aa,2));
        aa1 = bsxfun(@minus,aa,median(aaMod,1));
        Data1 = [Data1,aa1];
    end
    data = Data1;
end
data = (data - mean(data(:)))/sqrt(var(data(:)));
[nGenes nFeatures] = size(data);
%set parameters for proposal distribution for sampling of orders
uCT = unique(captureTimes);
n3 = floor(nFeatures/5);
n3a = floor(nFeatures/12);
n0 = floor(nFeatures/7);
delta = [1/4000,1/4000];
jj = 0.1;
kk = 0.1;
pp = [0,repmat(1/3,1,3)];
tau1 = sqrt(sum((diff(data,1,2).^2),1));
featureNames = cumsum([0 tau1])/sum(tau1);
%permute the data without capture times to obtain a random starting order
if permuteData == true
    xx = 1:nFeatures;
    yy = [];
    uniqueCaptureTimes = unique(captureTimes);
    lu = length(uniqueCaptureTimes);
    for i = 1:lu
        xy = xx(captureTimes == uniqueCaptureTimes(i));
        yy = [yy randsample(xy,length(xy))];
    end
    data = data(:,yy);
    tau1 = sqrt(sum((diff(data,1,2).^2),1));
    featureNames = cumsum([0 tau1])/sum(tau1);
     cellFileName = [strtok(fileName, '.'),'_Permutation', num2str(uniqueIdentifier)];
    CellFileName = [cellFileName '.csv'];
    csvwrite(CellFileName,yy); 
end

% For convenience, store the data in a structure
dataStruct.data         = data;
dataStruct.featureNames = featureNames;
%these are the original pseudotimes
[X, Y] = meshgrid(featureNames);
dataStruct.timeDiffs    = (-(X - Y).^2);
nGenes                  = size(data,1);
nFeatures               = length(featureNames);
dataStruct.nGenes       = nGenes;
dataStruct.nFeatures    = nFeatures;
dataStruct.order        = 1 : nFeatures;
dataStruct.proposedOrder= 1 : nFeatures;


%compute the probabilities for choosing moves
distanceMatrix1           = bsxfun(@minus,data,permute(data,[1 3 2]));
distanceMatrix            = sum(reshape(abs(distanceMatrix1),size(data,1),[]));
moveProb1                 = exp(-delta(1)* distanceMatrix.^2);
moveProb1                 = moveProb1./sum(moveProb1);
moveProb1                 = moveProb1.^jj;
dataStruct.moveProb1      = moveProb1./sum(moveProb1);
dataStruct.cumMoveProb1   = cumsum(dataStruct.moveProb1);
moveProb2                 = exp(-delta(2)* distanceMatrix.^2);
moveProb2                 = moveProb2./sum(moveProb2);
moveProb2                 = moveProb2.^kk;
dataStruct.moveProb2      = moveProb2./sum(moveProb2);
dataStruct.cumMoveProb2   = cumsum( dataStruct.moveProb2 );
dataStruct.pp            = cumsum(pp);

% Initialise the clustering partition (structOfClusters)
logL = normrnd(log(0.5),b);
a   = normrnd(log(sqrt(0.5)),0.01);
a1  = betarnd(4,1);
eps = 0.1*randn(1)+log(0.5);
squaredHypers   = [exp(2*logL);exp(a*2)*3+exp(eps);exp(2*a)*(1-a1);...
    1;exp(2*a)*a1];
[structOfClusters clusterIDs] = feval(fHandle, dataStruct, 'init',squaredHypers);
% Initialise the hyperparameters

logPriorOfLogHypers = [log(normpdf(logL,log((0.5)),b));log(normpdf(a,log(sqrt(0.5)),0.01));
    log(betapdf(a1,4,1));log(normpdf(eps,sqrt(0.5),0.1))];

a0 = gammaPrior(1); b0 = gammaPrior(2);  % Parameters for gamma prior
alpha = gamrnd(a0,1/b0); %Draw the initial alpha from the prior


% We will keep a track of the hyperparameter acceptances
nHyperProposals     = [0,0,0,0]; 
nHyperAcceptances   = nHyperProposals;
nOrderProposals     = zeros(1,4);
nOrderAcceptances   = zeros(1,4);
% Initialise the output file
outFile = [saveFileName '.csv'];
outFileOrder = [saveFileNameOrder '.csv'];
outFileLk    = [saveFileNameLk '.csv'];

fid = fopen(outFile, 'wt');
fopen(outFileOrder, 'wt');
fopen(outFileLk, 'wt');
fclose(fid);
  

%lines 172-339  were taken from the DPMSysBio code with minor changes
for sampleNumber = 1:nSamples

     if(mod(sampleNumber,clusterAllocationFrequency) == 0)
     
            if(verbose && mod(sampleNumber,10)==0)
                disp(['Sample number = ', num2str(sampleNumber)]);
            end
            % We iterate over the genes:
            for i = 1:dataStruct.nGenes
                % Which clusters are currently occupied?
                occupiedClusterIDs  = unique(clusterIDs);
                nOccupiedClusters   = length(occupiedClusterIDs);
                % Find the cluster in which the current gene resides
                clusterNumber   = clusterIDs(i);
                % Where does this cluster label appear in occupiedClusterIDs ?
                occupiedClusterIndex    = find(occupiedClusterIDs == clusterNumber);
                % Pick out the data for the current gene
                dataForCurrentGene      = data(i,dataStruct.order);
                % Pick out the occupied clusters
                currentClusters         = structOfClusters(occupiedClusterIDs);
                % For each occupied cluster, we will propose adding the gene to
                % each cluster (except for the cluster in which the current gene
                % currently resides, for which we will propose to remove the gene)
                proposedClusters        = currentClusters;
                % The following is a switch to indicate if, by removing the current
                % gene from its cluster, we have emptied the cluster
                emptiedClusterSwitch = false;
                for j = (1:nOccupiedClusters)
                    % We iterate through the clusters:
                    currentClusterLabel = occupiedClusterIDs(j);
                    currentProposedCluster = proposedClusters(j);
                    if( currentClusterLabel ~= clusterNumber)
                        % Then we need to add the gene to the cluster
                        nGenesInCluster   = currentProposedCluster.nGenes + 1;
                        dataCounts        = currentProposedCluster.dataCounts + dataForCurrentGene;
                        dataInCluster     = [currentProposedCluster.dataInCluster;dataForCurrentGene];
                        currentProposedClusterLogicalGeneIDs    = currentProposedCluster.logicalGeneIDs;
                        currentProposedClusterLogicalGeneIDs(i) = true;
                        currentProposedCluster.logicalGeneIDs   = currentProposedClusterLogicalGeneIDs;
                    else
                        % We need to remove the gene
                        nGenesInCluster   = currentProposedCluster.nGenes - 1;
                        dataCounts        = currentProposedCluster.dataCounts - dataForCurrentGene;
                        currentProposedClusterLogicalGeneIDs    = currentProposedCluster.logicalGeneIDs;
                        currentProposedClusterLogicalGeneIDs(i) = false;
                        currentProposedCluster.logicalGeneIDs   = currentProposedClusterLogicalGeneIDs;
                        dataInCluster =   data(currentProposedClusterLogicalGeneIDs,dataStruct.order);
                    end
                    % Store how many genes are left in the current gene's current
                    % cluster:
                    currentProposedCluster.nGenes            = nGenesInCluster;
                    if(nGenesInCluster > 0)%...then the cluster is not empty
                        % Update all of the required fields for the cluster
                        currentProposedCluster.dataCounts        = dataCounts;
                        currentProposedCluster.N                 = nGenesInCluster*nFeatures;
                        currentProposedCluster.dataInCluster     = dataInCluster;
                        if(length(currentProposedCluster.covarianceMatrixInverses) < nGenesInCluster || ...
                            isempty(currentProposedCluster.covarianceMatrixInverses(nGenesInCluster).determinant))
                            currentProposedCluster = feval(fHandle, currentProposedCluster, 'invert',squaredHypers); 
                        end
                        currentProposedCluster = feval(fHandle, currentProposedCluster, 'marginal',squaredHypers); 
                    else
                        % We must have emptied a cluster by removing the gene
                        emptiedClusterSwitch = true;
                    end
                    % Store the currentProposedCluster
                    proposedClusters(j) = currentProposedCluster;
                end
                if(emptiedClusterSwitch)
                   
                    % Then we just set the auxiliary component to be the component
                    % we just emptied
                    auxiliaryCluster = currentClusters(occupiedClusterIndex);
                    
                else                  
                    % We need to initialise the auxiliary component
                    auxiliaryCluster      = structOfClusters(end); %(an empty cluster)
                    %%we need the pseudotime order
                    auxiliaryCluster.timeDiffs = dataStruct.timeDiffs;
                    auxiliaryCluster      = feval(fHandle, auxiliaryCluster, 'initialiseAuxiliary',squaredHypers);
                    auxiliaryCluster.dataCounts        = dataForCurrentGene;
                    auxiliaryCluster.dataCountsProposedOrder        = dataForCurrentGene;
                    auxiliaryClusterLogicalGeneIDs    = false(1, nGenes);
                    auxiliaryClusterLogicalGeneIDs(i) = true;
                    auxiliaryCluster.logicalGeneIDs = auxiliaryClusterLogicalGeneIDs;
                    auxiliaryCluster.dataInCluster = dataForCurrentGene;
                    auxiliaryCluster.nGenes = 1;
                    auxiliaryCluster = feval(fHandle, auxiliaryCluster, 'marginal',squaredHypers);
                  
                end
                logMarginalLikelihoodsWithoutGene = [currentClusters.logMarginalLikelihood];
                logMarginalLikelihoodsWithGene    = [proposedClusters.logMarginalLikelihood];

                % We have to swap the entries for the current cluster (since we
                % removed rather than added the gene for this cluster)
                saved = logMarginalLikelihoodsWithoutGene(occupiedClusterIndex);
                logMarginalLikelihoodsWithoutGene(occupiedClusterIndex) = logMarginalLikelihoodsWithGene(occupiedClusterIndex);
                logMarginalLikelihoodsWithGene(occupiedClusterIndex)    = saved;

                % Calculate the log marginal likelihood ratios
                logMarginalLikelihoodRatios = logMarginalLikelihoodsWithGene - logMarginalLikelihoodsWithoutGene;
                marginalLikelihoodRatios    = exp(logMarginalLikelihoodRatios);

                geneCounts = [currentClusters.nGenes];
                % Again, we have to swap the entries for the current cluster
                geneCounts(occupiedClusterIndex) = proposedClusters(occupiedClusterIndex).nGenes;

                % Calculate the conditional posterior probabilities (unnormalised)
                existingClusterProbs  = geneCounts.*marginalLikelihoodRatios;
                newClusterProb        = alpha*exp(auxiliaryCluster.logMarginalLikelihood);

                % Normalise these probabilities
                allProbs   = [existingClusterProbs, newClusterProb];
                allProbs   = allProbs/sum(allProbs);

                % Determine the cluster to which the gene should be added:
                cumulProbs = cumsum(allProbs);

                selected = find(cumulProbs > rand, 1);
                % Update the clusters:
                if(selected~=occupiedClusterIndex)  % We only need to do anything if the gene has been put back in a different cluster
                    structOfClusters(clusterNumber) = proposedClusters(occupiedClusterIndex);
                    if(selected > nOccupiedClusters)% a new cluster is born...
                        % Get a new ID
                        available = setdiff((1:(nOccupiedClusters+1)), occupiedClusterIDs);
                        updatedClusterLabel = available(1);
                        if(updatedClusterLabel == (nOccupiedClusters+1))
                            % Extend the structre array by 1 element
                            structOfClusters(end+1) = structOfClusters(end);
                        end
                        
                        structOfClusters(updatedClusterLabel) = auxiliaryCluster;
                    else
                        % Find the appropriate existing ID
                        updatedClusterLabel = occupiedClusterIDs(selected);
                        structOfClusters(updatedClusterLabel) = proposedClusters(selected);

                    end
                    clusterIDs(i) = updatedClusterLabel;
                end
            end       
        k = length(unique(clusterIDs));
        eta = betarnd((alpha+1), nGenes);
        A = a0 + k - 1;
        B = nGenes*(b0 - log(eta));
        pi_eta = A/(A + B);
        if rand < pi_eta
        alpha = gamrnd(a0+k , 1/(b0 - log(eta)));
        else
        alpha = gamrnd(a0+k-1, 1/(b0 - log(eta)));
        end
     
        if(verbose && mod(sampleNumber,10)==0)
            disp(['alpha = ', num2str(alpha)]);
            disp(['nClusters =', num2str(length(unique(clusterIDs)))])
            disp(['Order acceptance rate = ', num2str(nOrderAcceptances./nOrderProposals)]);
        end
    end  
        
    if(mod(sampleNumber,thinningFreq) == 0)
        occupiedClusterIDs  = unique(clusterIDs);
        outputVector = [alpha, clusterIDs, squaredHypers'];
        outputOrder  = dataStruct.order;
        outputLk = sum([structOfClusters(occupiedClusterIDs).logMarginalLikelihood]/nGenes);
        dlmwrite(outFile, outputVector, '-append', 'delimiter', ',');
        dlmwrite(outFileOrder,outputOrder, '-append', 'delimiter',',');
        dlmwrite(outFileLk,outputLk, '-append', 'delimiter',',');
    end 
    %%% sample the hyperparameters
    if(mod(sampleNumber,hyperParameterSamplingFrequency) == 0)
        occupiedClusterIDs  = unique(clusterIDs);
        nOccupiedClusters   = length(occupiedClusterIDs);  
        %structure of proposed clusters
        proposedStructOfClusters = structOfClusters;
        emptyCovMatStruct = structOfClusters(end).covarianceMatrixInverses;
        % first sample L        
        proposedSquaredHypers   = squaredHypers;
        logLProp   = logL + randn*stepSize(1);
        currentLogPriorOfLogHypers  = logPriorOfLogHypers;
        proposedLogPriorOfLogHypers(1) = log(normpdf(logLProp,log(0.5),b));
        proposedSquaredHypers(1) = exp(2*logLProp); 
        sumCurrentMarginalLogLikelihoods  = 0.0;
        sumProposedMarginalLogLikelihoods = 0.0;
        nHyperProposals(1)  = nHyperProposals(1) + 1; 
        for i = 1:nOccupiedClusters
            clusterNumber    = occupiedClusterIDs(i);
            currentCluster   = structOfClusters(clusterNumber);
            currentCluster = feval(fHandle, currentCluster, 'marginal',squaredHypers);
            currentLogMarginalLikelihood    = currentCluster.logMarginalLikelihood;  
            sumCurrentMarginalLogLikelihoods = sumCurrentMarginalLogLikelihoods + ...
                currentLogMarginalLikelihood;
            %the cluster with the new GP parameters
            proposedCluster   = currentCluster;
            proposedCluster.covarianceMatrixInverses = emptyCovMatStruct;

            l2          = proposedSquaredHypers(1);
            sf2         = proposedSquaredHypers(2);
            se2         = proposedSquaredHypers(3);
            L2          = proposedSquaredHypers(4);
            SF2         = proposedSquaredHypers(5);

            timeDiffs   = proposedCluster.timeDiffs;                
            proposedCluster.covarianceMatrix        = sf2*exp(timeDiffs/(2*l2));
            proposedCluster.covarianceMatrix2       = SF2*exp(timeDiffs/(2*L2));
            proposedCluster.covarianceMatrix2Chol   = chol(proposedCluster.covarianceMatrix2 + se2*eye(proposedCluster.nFeatures));               

            proposedCluster = feval(fHandle, proposedCluster, 'invert',squaredHypers);
            proposedCluster = feval(fHandle, proposedCluster, 'marginal',squaredHypers);
            proposedLogMarginalLikelihood = proposedCluster.logMarginalLikelihood; 
            proposedStructOfClusters(clusterNumber) = proposedCluster;
            sumProposedMarginalLogLikelihoods = sumProposedMarginalLogLikelihoods + ...
                 proposedLogMarginalLikelihood; 

        end
        logRatio = sumProposedMarginalLogLikelihoods + ...
        sum(proposedLogPriorOfLogHypers) - ...
        sumCurrentMarginalLogLikelihoods - ...
        sum(currentLogPriorOfLogHypers); 
        %if we accept the new hyperparameter, we then update
        %the cluster structure
        if rand < exp(logRatio)
            structOfClusters = proposedStructOfClusters;
            logL = logLProp;
            squaredHypers = proposedSquaredHypers;
            logPriorOfLogHypers = proposedLogPriorOfLogHypers;
            nHyperAcceptances(1) = nHyperAcceptances(1) + 1;              
        end  
        %now sample a
        aProp                   = a + randn(1)*stepSize(2);                   
        currentLogPriorOfLogHypers      = logPriorOfLogHypers;
        proposedLogPriorOfLogHypers     = logPriorOfLogHypers;
        proposedLogPriorOfLogHypers(2)  = log(normpdf(aProp,log(sqrt(0.5)),0.01));
        proposedSquaredHypers       = squaredHypers; 
        proposedSquaredHypers(2)    = exp(aProp*2)*3+exp(eps);
        proposedSquaredHypers(3)    = exp(2*aProp)*(1-a1);
        proposedSquaredHypers(5)    = exp(2*(aProp))*a1;
        sumCurrentMarginalLogLikelihoods  = 0.0;
        sumProposedMarginalLogLikelihoods = 0.0;
        nHyperProposals(2)  = nHyperProposals(2) + 1; 
        for i = 1:nOccupiedClusters
            clusterNumber    = occupiedClusterIDs(i);
            currentCluster   = structOfClusters(clusterNumber);
            currentCluster = feval(fHandle, currentCluster, 'marginal',squaredHypers);
            currentLogMarginalLikelihood    = currentCluster.logMarginalLikelihood;  
            sumCurrentMarginalLogLikelihoods = sumCurrentMarginalLogLikelihoods + ...
                currentLogMarginalLikelihood;
            %the cluster with the new GP parameters
            proposedCluster   = currentCluster;
            proposedCluster.covarianceMatrixInverses = emptyCovMatStruct;

            l2          = proposedSquaredHypers(1);
            sf2         = proposedSquaredHypers(2);
            se2         = proposedSquaredHypers(3);
            L2          = proposedSquaredHypers(4);
            SF2         = proposedSquaredHypers(5);

            timeDiffs   = proposedCluster.timeDiffs;                
            proposedCluster.covarianceMatrix        = sf2*exp(timeDiffs/(2*l2));
            proposedCluster.covarianceMatrix2       = SF2*exp(timeDiffs/(2*L2));
            proposedCluster.covarianceMatrix2Chol   = chol(proposedCluster.covarianceMatrix2 + se2*eye(proposedCluster.nFeatures));               

            proposedCluster = feval(fHandle, proposedCluster, 'invert',squaredHypers);
            proposedCluster = feval(fHandle, proposedCluster, 'marginal',squaredHypers);
            proposedLogMarginalLikelihood = proposedCluster.logMarginalLikelihood; 
            proposedStructOfClusters(clusterNumber) = proposedCluster;
            sumProposedMarginalLogLikelihoods = sumProposedMarginalLogLikelihoods + ...
                 proposedLogMarginalLikelihood; 

        end
        logRatio = sumProposedMarginalLogLikelihoods + ...
        sum(proposedLogPriorOfLogHypers) - ...
        sumCurrentMarginalLogLikelihoods - ...
        sum(currentLogPriorOfLogHypers); 
        %if we accept the new hyperparameter, we then update
        %the cluster structure
        if rand < exp(logRatio)
            structOfClusters = proposedStructOfClusters;
            squaredHypers = proposedSquaredHypers;
            logPriorOfLogHypers = proposedLogPriorOfLogHypers;
            nHyperAcceptances(2) = nHyperAcceptances(2) + 1;   
            a = aProp;
        end  
        
        %now sample a1
        proposedSquaredHypers       = squaredHypers;
        a1Prop                      = a1 + randn(1)*stepSize(3); 
        if a1Prop < 1 && a1Prop >0 %else the posterior is 0    
            currentLogPriorOfLogHypers  = logPriorOfLogHypers;
            proposedLogPriorOfLogHypers = logPriorOfLogHypers;
            proposedLogPriorOfLogHypers(3) = log(betapdf(a1Prop,4,1));
            proposedSquaredHypers(3) = exp(2*a)*(1-a1Prop);
            proposedSquaredHypers(5) = exp(2*(a))*a1Prop; 
            sumCurrentMarginalLogLikelihoods  = 0.0;
            sumProposedMarginalLogLikelihoods = 0.0;
            nHyperProposals(3)  = nHyperProposals(3) + 1; 
            for i = 1:nOccupiedClusters
                clusterNumber    = occupiedClusterIDs(i);
                currentCluster   = structOfClusters(clusterNumber);
                currentCluster = feval(fHandle, currentCluster, 'marginal',squaredHypers);
                currentLogMarginalLikelihood    = currentCluster.logMarginalLikelihood;  
                sumCurrentMarginalLogLikelihoods = sumCurrentMarginalLogLikelihoods + ...
                    currentLogMarginalLikelihood;
                %the cluster with the new GP parameters
                proposedCluster   = currentCluster;
                proposedCluster.covarianceMatrixInverses = emptyCovMatStruct;

                l2          = proposedSquaredHypers(1);
                sf2         = proposedSquaredHypers(2);
                se2         = proposedSquaredHypers(3);
                L2          = proposedSquaredHypers(4);
                SF2         = proposedSquaredHypers(5);

                timeDiffs   = proposedCluster.timeDiffs;                
                proposedCluster.covarianceMatrix        = sf2*exp(timeDiffs/(2*l2));
                proposedCluster.covarianceMatrix2       = SF2*exp(timeDiffs/(2*L2));
                proposedCluster.covarianceMatrix2Chol   = chol(proposedCluster.covarianceMatrix2 + se2*eye(proposedCluster.nFeatures));               

                proposedCluster = feval(fHandle, proposedCluster, 'invert',squaredHypers);
                proposedCluster = feval(fHandle, proposedCluster, 'marginal',squaredHypers);
                proposedLogMarginalLikelihood = proposedCluster.logMarginalLikelihood; 
                proposedStructOfClusters(clusterNumber) = proposedCluster;
                sumProposedMarginalLogLikelihoods = sumProposedMarginalLogLikelihoods + ...
                     proposedLogMarginalLikelihood; 

            end
            logRatio = sumProposedMarginalLogLikelihoods + ...
            sum(proposedLogPriorOfLogHypers) - ...
            sumCurrentMarginalLogLikelihoods - ...
            sum(currentLogPriorOfLogHypers); 
            %if we accept the new hyperparameter, we then update
            %the cluster structure
            if rand < exp(logRatio)
                structOfClusters = proposedStructOfClusters;
                squaredHypers = proposedSquaredHypers;
                logPriorOfLogHypers = proposedLogPriorOfLogHypers;
                nHyperAcceptances(3) = nHyperAcceptances(3) + 1;  
                a1 = a1Prop;
            end  
        end
        %now sample eps
        epsProp                 = eps + randn(1)*stepSize(4);                   
        currentLogPriorOfLogHypers  = logPriorOfLogHypers;
        proposedLogPriorOfLogHypers  = logPriorOfLogHypers;
        proposedLogPriorOfLogHypers(4) = log(normpdf(epsProp,log(0.5),0.1));
        proposedSquaredHypers = squaredHypers; 
        proposedSquaredHypers(2) = exp(a*2)*3+exp(epsProp);
        sumCurrentMarginalLogLikelihoods  = 0.0;
        sumProposedMarginalLogLikelihoods = 0.0;
        nHyperProposals(4)  = nHyperProposals(4) + 1; 
        for i = 1:nOccupiedClusters
            clusterNumber    = occupiedClusterIDs(i);
            currentCluster   = structOfClusters(clusterNumber);
            currentCluster = feval(fHandle, currentCluster, 'marginal',squaredHypers);
            currentLogMarginalLikelihood    = currentCluster.logMarginalLikelihood;  
            sumCurrentMarginalLogLikelihoods = sumCurrentMarginalLogLikelihoods + ...
                currentLogMarginalLikelihood;
            proposedCluster   = currentCluster;
            proposedCluster.covarianceMatrixInverses = emptyCovMatStruct;

            l2          = proposedSquaredHypers(1);
            sf2         = proposedSquaredHypers(2);
            se2         = proposedSquaredHypers(3);
            L2          = proposedSquaredHypers(4);
            SF2         = proposedSquaredHypers(5);

            timeDiffs   = proposedCluster.timeDiffs;                
            proposedCluster.covarianceMatrix        = sf2*exp(timeDiffs/(2*l2));
            proposedCluster.covarianceMatrix2       = SF2*exp(timeDiffs/(2*L2));
            proposedCluster.covarianceMatrix2Chol   = chol(proposedCluster.covarianceMatrix2 + se2*eye(proposedCluster.nFeatures));               

            proposedCluster = feval(fHandle, proposedCluster, 'invert',squaredHypers);
            proposedCluster = feval(fHandle, proposedCluster, 'marginal',squaredHypers);
            proposedLogMarginalLikelihood = proposedCluster.logMarginalLikelihood; 
            proposedStructOfClusters(clusterNumber) = proposedCluster;
            sumProposedMarginalLogLikelihoods = sumProposedMarginalLogLikelihoods + ...
                 proposedLogMarginalLikelihood; 

        end
        logRatio = sumProposedMarginalLogLikelihoods + ...
        sum(proposedLogPriorOfLogHypers) - ...
        sumCurrentMarginalLogLikelihoods - ...
        sum(currentLogPriorOfLogHypers); 
        %if we accept the new hyperparameter, we then update
        %the cluster structure
        if rand < exp(logRatio)
            structOfClusters = proposedStructOfClusters;
            squaredHypers = proposedSquaredHypers;
            logPriorOfLogHypers = proposedLogPriorOfLogHypers;
            nHyperAcceptances(4) = nHyperAcceptances(4) + 1;   
            eps = epsProp;
        end  
        
        
        if(verbose && mod(sampleNumber,5*hyperParameterSamplingFrequency)==0)
            disp('Sampling hyperparameters...')
            disp(['Hypers acceptance rate = ', num2str(nHyperAcceptances./nHyperProposals)]);         
        end
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %now we sample the orders
   [dataStruct.proposedOrder dataStruct.proposedPseudoTimes move] = ProposingPseudoTimes(dataStruct,false,n0,n3,n3a);
    nOrderProposals(move+1) = nOrderProposals(move+1) + 1;
    occupiedClusterIDs  = unique(clusterIDs);
    nOccupiedClusters   = length(occupiedClusterIDs);
    %we order the data according to the proposed order
    reorderedData       = dataStruct.data(:,dataStruct.proposedOrder);
    
    LogMarginalLikelihoods              = zeros([1,nOccupiedClusters]);
    LogMarginalLikelihoodsProposedOrder = zeros([1,nOccupiedClusters]);
    [X, Y] = meshgrid(dataStruct.proposedPseudoTimes);
    timeDiffsProposedOrder   = (-(X - Y).^2);
    l2          = squaredHypers(1);
    sf2         = squaredHypers(2);
    L2          = squaredHypers(4);
    SF2         = squaredHypers(5);
    for i = 1:nOccupiedClusters
        clusterNumber                           = occupiedClusterIDs(i);
        currentCluster                          = structOfClusters(clusterNumber);
        currentCluster.timeDiffsProposedOrder   =  timeDiffsProposedOrder;
        currentCluster.dataInClusterProposedOrder         = ...
        reorderedData(currentCluster.logicalGeneIDs,:);        
        currentCluster.dataCountsProposedOrder  = sum(currentCluster.dataInClusterProposedOrder,1); 
        %we compute the covariance matrices corresponding to the proposed
        %order
        currentCluster.covarianceMatrixProposedOrder    = sf2*exp(currentCluster.timeDiffsProposedOrder/(2*l2));
        currentCluster.covarianceMatrix2ProposedOrder   = SF2*exp(currentCluster.timeDiffsProposedOrder/(2*L2));
        %we compute the required matrix decompositions and likelihoods for
        %the proposed order
        currentCluster = feval(fHandle, currentCluster, 'invertSamplingOrder',squaredHypers);        
        currentCluster = feval(fHandle, currentCluster, 'sampleOrders',squaredHypers);
        LogMarginalLikelihoods(i)               = currentCluster.logMarginalLikelihood;
        LogMarginalLikelihoodsProposedOrder(i)  = currentCluster.logMarginalLikelihoodProposedOrder;
        structOfClusters(clusterNumber) = currentCluster;
        
    end
    logRatio = sum(LogMarginalLikelihoodsProposedOrder) - sum(LogMarginalLikelihoods);
    %we accept/reject the proposed order in a MH step
    %if the proposed order is accepted we update the data struct of the
    %clusters
    if rand < exp(logRatio)
        dataStruct.order            = dataStruct.proposedOrder;
        dataStruct.featureNames     = dataStruct.proposedPseudoTimes;
        nOrderAcceptances(move+1)   = nOrderAcceptances(move+1) + 1;
        %se2 = min(1-sf2,0.1);
        se2 = squaredHypers(3);
        for i = 1:nOccupiedClusters
            clusterNumber    = occupiedClusterIDs(i);
            currentCluster   = structOfClusters(clusterNumber); 
            nGenesInCluster    = currentCluster.nGenes;
            currentCluster.dataCounts           = currentCluster.dataCountsProposedOrder;
            currentCluster.dataInCluster = currentCluster.dataInClusterProposedOrder;          
            currentCluster.covarianceMatrix = currentCluster.covarianceMatrixProposedOrder;
            currentCluster.covarianceMatrix2 = currentCluster.covarianceMatrix2ProposedOrder;
            currentCluster.covarianceMatrix2Chol = chol(currentCluster.covarianceMatrix2 + se2*eye(currentCluster.nFeatures));
            currentCluster.covarianceMatrixInverses(1,nGenes) =...
            struct('invertedCovarianceMatrix', [], 'determinant', []);
            currentCluster.covarianceMatrixInverses(nGenesInCluster).determinant=currentCluster.covarianceMatrixInversesProposedOrder(nGenesInCluster).determinant;   
            currentCluster.covarianceMatrixInverses(nGenesInCluster).invertedCovarianceMatrix = ...
            currentCluster.covarianceMatrixInversesProposedOrder(nGenesInCluster).invertedCovarianceMatrix; 
            currentCluster.SigmaInv = currentCluster.SigmaInvProposedOrder;  
            currentCluster.timeDiffs = currentCluster.timeDiffsProposedOrder;
            if i == 1
                dataStruct.timeDiffs = currentCluster.timeDiffs;
            end
            currentCluster.covarianceMatrix2Chol = chol(currentCluster.covarianceMatrix2 + squaredHypers(3)*eye(currentCluster.nFeatures));
            currentCluster.logMarginalLikelihood = currentCluster.logMarginalLikelihoodProposedOrder;
             currentCluster = feval(fHandle, currentCluster, 'marginal',squaredHypers);
             structOfClusters(clusterNumber) = currentCluster;  
        end
       
    end
end   


end
