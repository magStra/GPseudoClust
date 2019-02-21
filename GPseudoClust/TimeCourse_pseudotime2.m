function [output1 output2] = TimeCourse_pseudotime2(input, mode,squaredHypers)
switch mode
    case 'init'
        %initialises a data struct of clusters to store the current state
        %of the chain
        data              = input.data;
        featureNames      = input.featureNames;
        nGenes            = input.nGenes;
        nFeatures         = input.nFeatures;
        [X, Y] = meshgrid(featureNames);
        timeDiffs         = (-(X - Y).^2);
        nStartingClusters   = nGenes;
        % clusterIDs          = random('unid', nStartingClusters, 1, nGenes);
        clusterIDs = 1:nGenes; 
        uniqueIDs         = unique(clusterIDs);
        sparseMatrix      = zeros(nGenes,nFeatures);
       
        % Define the cluster structure
        clusterStruct(1,nGenes+2) = struct(...
            'nGenes', [], ...
            'nFeatures', [], ...
            'timeDiffs', [],...
            'timeDiffsProposedOrder',[],...
            'N', [], ...
            'dataCounts', [], ...
            'dataCountsProposedOrder',[],...
            'logicalGeneIDs', [], ...
            'logMarginalLikelihood',[], ...
            'logMarginalLikelihoodProposedOrder',[],...
            'covarianceMatrix', [], ...
            'covarianceMatrixProposedOrder',[],...
            'covarianceMatrixInverses', [], ...
            'covarianceMatrixInversesProposedOrder',[],...
            'covarianceMatrix2',[],...
            'covarianceMatrix2Chol',[],...
            'SigmaInv',[],...
            'SigmaInvProposedOrder',[],...
            'covarianceMatrix2ProposedOrder',[],...
            'dataInCluster',[],...
            'dataInClusterProposedOrder',[]);
        
        [clusterStruct.nFeatures       ] = deal(nFeatures);
        [clusterStruct.timeDiffs] = deal(timeDiffs);
     
        for i = uniqueIDs
            logicalIndices                   = clusterIDs == i;
            indices                          = find(logicalIndices);
            nGenesInCluster                  = length(indices);
            dataInCluster                    = sparseMatrix;
            dataInCluster(indices,:)         = data(logicalIndices,:);
            currentCluster                   = clusterStruct(i);
            dataInClusterInd                 = dataInCluster(indices,:);
            currentCluster.dataInCluster     = dataInClusterInd;
            currentCluster.dataInClusterProposedOrder     = dataInClusterInd;
            currentCluster.logicalGeneIDs    = logicalIndices;
            currentCluster.dataCounts        = sum(dataInCluster,1);
            currentCluster.dataCountsProposedOrder =  currentCluster.dataCounts;
            currentCluster.nGenes            = nGenesInCluster;
            currentCluster.N                 = nFeatures*nGenesInCluster;
                    
            l2  = squaredHypers(1);
            sf2 = squaredHypers(2);
            se2 = squaredHypers(3);
            L2  = squaredHypers(4);
            SF2 = squaredHypers(5);
           
            currentCluster.covarianceMatrix = sf2*exp(timeDiffs/(2*l2));%squared
            %exponential covariance matrix
            currentCluster.covarianceMatrixProposedOrder = currentCluster.covarianceMatrix;
            currentCluster.covarianceMatrix2 = SF2*exp(timeDiffs/(2*L2));
            currentCluster.covarianceMatrix2Chol = chol(currentCluster.covarianceMatrix2 + se2*eye(currentCluster.nFeatures));
            currentCluster.covarianceMatrix2ProposedOrder = currentCluster.covarianceMatrix2;

            currentCluster.covarianceMatrixInverses =...
                struct('invertedCovarianceMatrix', [], 'determinant', []);
            currentCluster.covarianceMatrixInversesProposedOrder =...
                struct('invertedCovarianceMatrixProposedOrder', [], 'determinantProposedOrder', []);
            
            currentCluster      = TimeCourse_pseudotime2(currentCluster, 'invert',...
                squaredHypers);
            currentCluster      = TimeCourse_pseudotime2(currentCluster, 'invertSamplingOrder',...
                squaredHypers);
            currentCluster = TimeCourse_pseudotime2(currentCluster, 'marginal',...
                squaredHypers);
            currentCluster = TimeCourse_pseudotime2(currentCluster, 'sampleOrders',...
                squaredHypers);
            clusterStruct(i) = currentCluster;
        end
        currentCluster = clusterStruct(nStartingClusters+1);
        currentCluster.covarianceMatrixInverses(1,nGenes) =...
            struct('invertedCovarianceMatrix', [], 'determinant', []);
        currentCluster.covarianceMatrixInversesProposedOrder(1,nGenes) =...
            struct('invertedCovarianceMatrixProposedOrder', [],'determinantProposedOrder', []);
        clusterStruct(nStartingClusters+1) = currentCluster;
        clusterStruct(end) = [];      
        output1 = clusterStruct;
        output2 = clusterIDs;
    
        
    case 'marginal'
        nGenesInCluster    = input.nGenes;
        N                  = input.N;
        dataCounter        = input.dataCounts;
        invertedMatrix     = input.covarianceMatrixInverses(nGenesInCluster).invertedCovarianceMatrix;
        logDetK            = input.covarianceMatrixInverses(nGenesInCluster).determinant;
        yColSum            = dataCounter';
        y2                 = input.dataInCluster;
        invertedSigma      = input.SigmaInv;
        abc = 0;
        for i = 1:nGenesInCluster
            abc = abc + y2(i,:)*invertedSigma*y2(i,:)';
            if isnan(abc)
                'NA'
            end
        end
        input.logMarginalLikelihood = -0.5*(N*log(2*pi)+logDetK+(yColSum'*invertedMatrix*yColSum + abc));
        output1  = input;
        output2  = [];
    case 'sampleOrders'
        % computes the contribution to the log marginal likelihood of one
        % cluster of genes after a new order has been proposed
        nGenesInCluster    = input.nGenes;
        N                  = input.N;
        invertedMatrix     = input.covarianceMatrixInversesProposedOrder(nGenesInCluster).invertedCovarianceMatrix;
        logDetK            = input.covarianceMatrixInversesProposedOrder(nGenesInCluster).determinant;
        yColSum            = input.dataCountsProposedOrder';
        y2                 = input.dataInClusterProposedOrder;
        invertedSigma      = input.SigmaInvProposedOrder;
        abc = 0;
        for i = 1:nGenesInCluster
            abc = abc + y2(i,:)*invertedSigma*y2(i,:)';
        end
        input.logMarginalLikelihoodProposedOrder = -0.5*(N*log(2*pi)+logDetK+(yColSum'*invertedMatrix*yColSum + abc));
        output1  = input;
        output2  = [];
    case 'invertSamplingOrder'
        %computes the matrix decompositions required for the likelihood
        %computations for proposed new orders
        nGenesInCluster    = input.nGenes;
        sf2                = squaredHypers(2);
        se2 = squaredHypers(3);
        nTimes             = input.nFeatures;       
        K                  = input.covarianceMatrixProposedOrder;
        Sigma              = input.covarianceMatrix2ProposedOrder;
        SigmaSE2           = Sigma + se2*eye(nTimes); 
        SigmaChol          = chol(SigmaSE2);
        opts.LT     = false;
        opts.UT     = true;
        X           = linsolve(SigmaChol,eye(nTimes),opts);
        SigmaInv    = X*X';
        try
            Chol2       = chol(nGenesInCluster*K + SigmaSE2);
        
        opts.LT     = true;
        opts.UT     = false;
       
        Y = linsolve(Chol2',-K*SigmaInv,opts);
        opts.LT = false;
        opts.UT = true;
        invertedK          = linsolve(Chol2,Y,opts);
        logDetK            =  2*(sum(log(diag(Chol2)))) + (nGenesInCluster-1)*2*sum(log(diag(SigmaChol)));
       
        if(~isreal(logDetK))  %We should not ever enter here, but just in case
            disp('Sampling orders: Numerical error - covariance matrix may not be positive definite')
        end
        
        input.covarianceMatrixInversesProposedOrder(nGenesInCluster).invertedCovarianceMatrix...
            = invertedK;
        input.covarianceMatrixInversesProposedOrder(nGenesInCluster).determinant = logDetK;
        input.SigmaInvProposedOrder = SigmaInv;
        catch
            'cholError'
        end
        output1 = input;
        output2  = [];
    
    case 'invert'
        %computes the matrix decompositions required for the likelihood
        %computations
        nGenesInCluster    = input.nGenes;
        sf2                = squaredHypers(2);
        se2 = squaredHypers(3);
        nTimes             = input.nFeatures;
        K                  = input.covarianceMatrix;
        Sigma              = input.covarianceMatrix2;
        SigmaSE2           = Sigma + se2*eye(nTimes); 
        SigmaChol          = input.covarianceMatrix2Chol;
        opts.LT     = false;
        opts.UT     = true;
        X           = linsolve(SigmaChol,eye(nTimes),opts);
        SigmaInv    = X*X';
        opts.LT     = true;
        opts.UT     = false;
       
        try
            Chol2       = chol(nGenesInCluster*K + SigmaSE2);
       
        Y = linsolve(Chol2',-K*SigmaInv,opts);
        opts.LT = false;
        opts.UT = true;
        invertedK          = linsolve(Chol2,Y,opts);
        logDetK            =  2*(sum(log(diag(Chol2)))) + (nGenesInCluster-1)*2*sum(log(diag(SigmaChol)));
     
        if(~isreal(logDetK))  %We should not ever enter here, but just in case
            disp('Numerical error - covariance matrix may not be positive definite')
        end
        input.SigmaInv = SigmaInv;
        input.covarianceMatrixInverses(nGenesInCluster).invertedCovarianceMatrix...
            = invertedK;
        input.covarianceMatrixInverses(nGenesInCluster).determinant = logDetK;
        input.SigmaInv = SigmaInv;
        catch
            'cholError'
        end
        output1 = input;
        output2  = [];
        
    case 'initialiseAuxiliary'
        %initialise an auxiliary cluster (Neal 2000, algorithm 8)
        output1 = input;
        nGenesInCluster    = 1;
        timeDiffs          = output1.timeDiffs;
        output1.nGenes     = nGenesInCluster;
        output1.N          = input.nFeatures;
        l2                    = squaredHypers(1);
        sf2                   = squaredHypers(2);
       se2 = squaredHypers(3);
        L2                    = squaredHypers(4);
        SF2                   = squaredHypers(5);
        output1.covarianceMatrix = sf2*exp(timeDiffs/(2*l2));
        output1.covarianceMatrix2 = SF2*exp(timeDiffs/(2*L2));
        output1.covarianceMatrix2Chol = chol(output1.covarianceMatrix2 + se2*eye(output1.nFeatures));
        output1      = TimeCourse_pseudotime2(output1, 'invert',squaredHypers);
        output2 = [];
end






