function [] = simClustLarge(uniqueIdentifier)
    addpath(genpath('../'));
    set(0,'defaultFigureUnits','centimeters');
    inputSeed = sum(100*clock) + sum(100*double(uniqueIdentifier)); 
    sim2Large(sprintf('simLarge%d',uniqueIdentifier),inputSeed);
    fHandle          =  @GPseudoClust2;  
    fileName         = sprintf('simLarge%d.csv',uniqueIdentifier);
    nSamples         = 5000;  
    verbose          = false; 
    inputSeed        = NaN;    
    permuteData      = true;
    b = 0.01;
    adjustForCellSize = false;
    A = csvread(fileName);
    [nGenes nCells] = size(A);
    a = nCells/3;
    PSMs = zeros(nGenes,nGenes,96);
    
    for k = 1:3
        captureTimes     = [ones([1,10*k]),repmat(2,[1,10*k]),repmat(3,[1,10*k])];
        parpool(12);
        tic
        for jj = 1:96
         ab = [randsample(1:a,10*k),randsample((a+1):(2*a),10*k),randsample((2*a+1):(3*a),10*k)];
         feval(fHandle, fileName, jj, nSamples, verbose,inputSeed,...
                        permuteData,captureTimes,b,adjustForCellSize,ab);
         PSMs(:,:,jj) = psm(dlmread(sprintf('simLarge%d_Results_Chain%d.csv',uniqueIdentifier,jj),',',[2449 1 4999 nGenes]));
        end
        toc
        pool = gcp('nocreate');
        delete(pool);
        save(sprintf('PSMsLarge%d_%d.mat',10*k,uniqueIdentifier),'PSMs');
    end 
end

%This was run with identifiers 1,2,3,..., 10, [code: simClustLarge(1), ...
%simClustLarge(10), which generated 10 simulated data sets and ran GPSseudoClust 
%with subsampled chains with 10, 20 and 30 cells per capture time for each
%of the simulated data sets. 

% %code to plot the generated data sets with smaller dots
% set(0,'defaultFigureUnits','centimeters');
% set(0,'DefaultAxesFontSize',20);
% nCells = 9000;
% tau = (0.5:(nCells-0.5))/nCells;
% for k = 1:10
%     simData = csvread(sprintf('simLarge%d.csv',k));
%     load(sprintf('simLarge%d.mat',k));
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
%      print(sprintf('simLarge%d.pdf',k),'-dpdf');
% end

%computing times - 12 cores for 96 subsampled chains

% Elapsed time is 3219.912705 seconds. (Wilkes, slower than BSU-CPU cluster)
% Elapsed time is 5771.455169 seconds.
% Elapsed time is 17271.756433 seconds.
% 
% Elapsed time is 2900.760924 seconds.(Wilkes)
% Elapsed time is 6295.033491 seconds.
% Elapsed time is 13145.550876 seconds.
% 
% Elapsed time is 2199.192744 seconds.(Wilkes)
% Elapsed time is 4673.942223 seconds.
% Elapsed time is 10371.339753 seconds.
% 
% Elapsed time is 4896.852603 seconds. (Wilkes)
% Elapsed time is 9659.345137 seconds.
% Elapsed time is 17737.512525 seconds.
% 
% Elapsed time is 2894.683964 seconds. (Wilkes)
% Elapsed time is 4354.596135 seconds.
% Elapsed time is 10992.540598 seconds.
% 
% avWilkes = mean([3219.912705,2900.760924,2199.192744,4896.852603,2894.683964]);
% avWilkes20 = mean([5771.455169,6295.033491,4673.942223,9659.345137,4354.596135]);
% avWilkes30 = mean([17271.756433,13145.550876,10371.339753,17737.512525,10992.540598]);
% 
% Elapsed time is 2683.971107 seconds.(BSU-CPU)
% Elapsed time is 4344.899608 seconds.
% Elapsed time is 10507.376108 seconds.
% 
% Elapsed time is 2671.261665 seconds.(BSU-CPU)
% Elapsed time is 3936.384937 seconds.
% Elapsed time is 7953.896243 seconds.
% 
% Elapsed time is 3582.380049 seconds.(BSU-CPU)
% Elapsed time is 4835.900080 seconds.
% Elapsed time is 11500.911111 seconds.
% 
% Elapsed time is 2401.934442 seconds. (BSU-CPU)
% Elapsed time is 3689.505829 seconds.
% Elapsed time is 8571.786738 seconds.
% 
% Elapsed time is 1705.536820 seconds.(BSU-CPU)
% Elapsed time is 3048.714884 seconds.
% Elapsed time is 5948.961344 seconds.

 %median compute times for BSU-CPU: 10 cells per capture time: 2671 sec
 %20 cells: 3936 sec, 30 cells: 8572 sec

 
 %code to read in samples of concentration parameters from output files
% alphas = [];
% for k = 1:10
%     for j = 1:96
%         File = sprintf('simLarge%d_Results_Chain%d.csv',k,j);
%         xx = dlmread(File,',',[0 0 2500-1 0]);
%         alphas = [alphas xx];
%     end
%     csvwrite(sprintf('alphasSimLarge%d.csv',k),alphas);
% end



