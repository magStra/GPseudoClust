% GPseudoClust applied to second simulated data set (less clearly separated
% clusters) with additional dropout noise
%100 simulations
%for 60% of the subsamples, all the genes are affected by dropout, for 40%
%only 50% of the genes are affected
addpath(genpath(pwd));
%add path to subaxis function
%https://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot
addpath(genpath('~/subaxis'));
for j = 1:100
    j
    A  = csvread('simDataClust2.csv');
    %select the 60% of cells with all genes affected by dropout, and the
    %40% with only 50% affected
    set1 = [randsample(1:20,12),randsample(21:40,12),randsample(41:60,12)];
    set2 = setdiff(1:60,set1);
    for k = 1:52
        nn = randsample(1:20,1);%up to 50% of the selected 40 cells are affected 
        xx = randsample(set1,nn);
        A(k,xx) = 0;
    end
    %for the other 40% of cells only 50% of the genes are affected
    for k = set2
        %affected with a uniform [0,0.5] probability
        prob = 0.5*rand();
        if rand() < prob
            yy = randsample(1:52,26);
            A(yy,k) = 0;
        end
    end
    csvwrite('simDataClustMixed.csv',A);
    captureTimes     = [ones([1,5]),repmat(2,[1,5]),repmat(3,[1,5])];
    %3 simulated capture times
    fHandle          =  @GPseudoClust2;  
    fileName         = 'simDataClustMixed.csv';
    uniqueIdentifiers = 1:24;
    nSamples         = 1000;  
    verbose          = false; 
    inputSeed        = NaN;    
    permuteData      = true;
    b                = 0.01;
    adjustForCellSize = false;
    tic
    parfor jk = 1:24
        %14 of the subsets are from set1, 10 from set 2
        if jk < 15
            subS = [randsample(intersect(1:20,set1),5),...
            randsample(intersect(21:40,set1),5),randsample(intersect(41:60,set1),5)];
        else
            subS = [randsample(intersect(1:20,set2),5),...
            randsample(intersect(21:40,set2),5),randsample(intersect(41:60,set2),5)];
        end
        feval(fHandle, fileName, uniqueIdentifiers(jk), nSamples,verbose, inputSeed,...
        permuteData,captureTimes,b,adjustForCellSize,subS);
    end
    toc
   %posterior similarity matrix
   PSMs = zeros(52,52,24);
   for jk = 1:24
       B = dlmread(sprintf('simDataClustMixed_Results_Chain%d.csv',jk),',',[500 1 999 52]);
       PSMs(:,:,jk) = psm(B);
   end
   save(sprintf('PSM_simDropout_Mixed_%d.mat',j),'PSMs');
   delete *simDataClustMixed_*
end

clear all;
%processing the PSMs
%using localised kernel k-means
addpath(genpath(pwd));
sumClustSimDropout_llkm = zeros(52,100);
for j = 1:100
    j
    K = load(sprintf('~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/dropout/PSM_simDropout_Mixed_%d.mat',j)); 
    K = K.PSMs;
    clusterSolution = computeSummaryPSM_lmkk(K,2:10);
    csvwrite(sprintf('~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/dropout/PSM_simDropoutMixed_lmkk_%d.csv',j),...
        clusterSolution.PSM);
    sumClustSimDropout_llkm(:,j) = clusterSolution.clustering;
end
csvwrite('~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/dropout/sumClust_simDropoutMixed_lmkk.csv',...
     sumClustSimDropout_llkm);

 clear all;
% consensus method
for j = 1:100
    j
    K = load(sprintf('~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/dropout/PSM_simDropout_Mixed_%d.mat',j)); 
    K = K.PSMs;
    PSM = mean(K,3);
    csvwrite(sprintf('~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/dropout/PSM_simDropoutMixed_consensus_%d.csv',j),...
        PSM);
end

figure()
a = randsample(1:100,12,false);
for j = 1:12
 
    K = csvread(sprintf('~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/dropout/PSM_simDropoutMixed_lmkk_%d.csv',a(j)));
     ax(j) =  subaxis(4,3,j, 'SpacingHoriz', 0.03, 'SpacingVert',0.03,'Padding', 0.01, 'MarginRight', 0.1,...
         'MarginLeft',0.05,'MarginTop',0.01,'MarginBottom',0.05);
    imagesc(K);
    colormap(hot);
     set(gca,'YDir','normal')
end
s2Pos = get(ax(12),'position');
s1Pos = get(ax(3),'position');
col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.01  s2Pos(2)  0.05  s1Pos(2)+s1Pos(3)*0.56])
col.FontSize = 12;
set(gcf, 'PaperUnits', 'centimeters');
 x_width=21.7 ;y_width=23;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('PSMS_dropoutMixed_lmkk','-dpdf','-r300')

figure()
for j = 1:12
 
     K = csvread(sprintf('~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/dropout/PSM_simDropoutMixed_consensus_%d.csv',a(j)));      

     ax(j) =  subaxis(4,3,j, 'SpacingHoriz', 0.03, 'SpacingVert',0.03,'Padding', 0.01, 'MarginRight', 0.1,...
         'MarginLeft',0.05,'MarginTop',0.01,'MarginBottom',0.05);
    imagesc(K);
    colormap(hot);
     set(gca,'YDir','normal')
end
s2Pos = get(ax(12),'position');
s1Pos = get(ax(3),'position');
col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.01  s2Pos(2)  0.05  s1Pos(2)+s1Pos(3)*0.56])
col.FontSize = 12;
set(gcf, 'PaperUnits', 'centimeters');
 x_width=21.7 ;y_width=23;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('PSMS_dropoutMixed_consensus','-dpdf','-r300')


