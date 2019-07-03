function [] = simClustLargeMoreNoise(uniqueIdentifier)
    addpath(genpath('../'));
    set(0,'defaultFigureUnits','centimeters');
    inputSeed = sum(100*clock) + sum(100*double(uniqueIdentifier)); 
    sim2LargeMoreNoise(sprintf('simLargeMoreNoise%d',uniqueIdentifier),inputSeed);
    fHandle          =  @GPseudoClust2;  
    fileName         = sprintf('simLargeMoreNoise%d.csv',uniqueIdentifier);
    nSamples         = 5000;  
    verbose          = false; 
    inputSeed        = NaN;    
    permuteData      = true;
    b = 0.01;
    adjustForCellSize = false;
    A = csvread(fileName);
    [nGenes nCells] = size(A);
    a = nCells/3;
    PSMs = zeros(nGenes,nGenes,48);
    
    for k = 1:3
        captureTimes     = [ones([1,10*k]),repmat(2,[1,10*k]),repmat(3,[1,10*k])];
        parpool(12);
        tic
        parfor jj = 1:96
         ab = [randsample(1:a,10*k),randsample((a+1):(2*a),10*k),randsample((2*a+1):(3*a),10*k)];
         feval(fHandle, fileName, jj, nSamples, verbose,inputSeed,...
                        permuteData,captureTimes,b,adjustForCellSize,ab);
         PSMs(:,:,jj) = psm(dlmread(sprintf('simLargeMoreNoise%d_Results_Chain%d.csv',uniqueIdentifier,jj),',',[2449 1 4999 nGenes]));
        end
        toc
        pool = gcp('nocreate');
        delete(pool);
        save(sprintf('PSMsLargeMoreNoise%d_%d.mat',10*k,uniqueIdentifier),'PSMs');
    end 
end

%This was run with identifiers 1,2, which generated 2 simulated data sets and ran GPSseudoClust 
%with subsampled chains with 10, 20 and 30 cells per capture time for each
%of the simulated data sets. 

%on BSU-Wilkesfor 96 subsampled chains (simulated data sets 1 and 2)
%Elapsed time is 2008.662524 seconds.(10 cells)
%Elapsed time is 5270.617253 seconds. (20 cells)
%Elapsed time is 7791.506486 seconds. (30 cells)

%With 48 subsampled chains again on BSU-Wilkes
%Elapsed time is 1204.427843 seconds. (10 cells)
%Elapsed time is 1618.268104 seconds.
%Elapsed time is 1330.244122 seconds.
%Elapsed time is 1948.924937 seconds. (20 cells)
%Elapsed time is 2732.914316 seconds.
%Elapsed time is 1905.230293 seconds.
%Elapsed time is 9209.582797 seconds. (30 cells)
%Elapsed time is 7152.547214 seconds.
%Elapsed time is 3426.797461 seconds.



% %code to plot the generated data sets with smaller dots
% set(0,'defaultFigureUnits','centimeters');
% set(0,'DefaultAxesFontSize',20);
% nCells = 9000;
% tau = (0.5:(nCells-0.5))/nCells;
% for k = 1:2
%     simData = csvread(sprintf('simLargeMoreNoise%d.csv',k));
%     load(sprintf('simLargeMoreNoise%d.mat',k));
%     colors = ['g','m','b','r','k'];
%     set(0,'defaultfigurecolor',[1 1 1])
%     set(0,'defaultaxesfontsize',12);
%     figure()
%     for j = 1:nGenes
%         plot(tau,simData(j,:),'o','Color',colors(z(j)),'MarkerSize',0.3);
%         hold on;
%     end
%      set(gcf, 'PaperPosition', [0 0 8.8 6]);
%      set(gcf,'PaperSize',[8.8 6]);
%      print(sprintf('simLargeMoreNoise%d.pdf',k),'-dpdf');
% end

%code to read in samples of concentration parameters from output files
% alphas = [];
% for k = 1:5
%     for j = 1:48
%         File = sprintf('simLargeMoreNoise%d_Results_Chain%d.csv',k,j);
%         xx = dlmread(File,',',[0 0 2500-1 0]);
%         alphas = [alphas xx];
%     end
%     csvwrite(sprintf('alphasSimLargeMoreNoise%d.csv',k),alphas);
% end
