% GPseudoClust applied to second simulated data set (less clearly separated
% clusters) with additional dropout noise (higher levels of dropout)
%100 simulations
addpath(genpath('../'));
%add path to subaxis function
%https://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot
addpath(genpath('~/subaxis'));
for j = 1:100
    j
    A  = csvread('simDataClust2.csv');
    %add dropout noise
    for k = 1:52
        nn = randsample(1:30,1);
        xx = randsample(1:60,nn);
        A(k,xx) = 0;
    end
    csvwrite('simDataClust2a.csv',A);
    captureTimes     = [ones([1,10]),repmat(2,[1,10]),repmat(3,[1,10])];
    %3 simulated capture times
    fHandle          =  @GPseudoClust2;  
    fileName         = 'simDataClust2a.csv';
    uniqueIdentifiers = 1:24;
    nSamples         = 1000;  
    verbose          = false; 
    inputSeed        = NaN;    
    permuteData      = true;
    b                = 0.01;
    adjustForCellSize = false;
    tic
    parfor jk = 1:24
        subS = [randsample(1:20,10),randsample(21:40,10),randsample(41:60,10)];
        feval(fHandle, fileName, uniqueIdentifiers(jk), nSamples,verbose, inputSeed,...
        permuteData,captureTimes,b,adjustForCellSize,subS);
    end
    toc
   %posterior similarity matrix
   PSMs = zeros(52,52,24);
   for jk = 1:24
       B = dlmread(sprintf('simDataClust2a_Results_Chain%d.csv',jk),',',[500 1 999 52]);
       PSMs(:,:,jk) = psm(B);
   end
   save(sprintf('PSM_simDropout30_Subsample_%d.mat',j),'PSMs');
   delete *simDataClust2a_*
end

%processing the PSMs
%using localised kernel k-means
addpath(genpath(pwd));
sumClustSimDropout30_llkm = zeros(52,100);
for j = 1:100
    j
    K = load(sprintf('PSM_simDropout30_Subsample_%d.mat',j)); 
    K = K.PSMs;
    clusterSolution = computeSummaryPSM_lmkk(K,2:10);
    csvwrite(sprintf('PSM_simDropout30_lmkk_%d.csv',j),...
        clusterSolution.PSM);
  
    sumClustSimDropout30_llkm(:,j) = clusterSolution.clustering;

end
csvwrite('sumClust_simDropout30_lmkk.csv',sumClustSimDropout30_llkm);
 
figure()
for j = 1:12
 
     K = csvread(sprintf('PSM_simDropout_consensus30_%d.csv',a(j)));      

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
print('PSMS_dropout30_consensus','-dpdf','-r300')

figure()
for j = 1:12
    load(sprintf('PSM_simDropout30_DPM_%d.mat',a(j)));
     ax(j) =  subaxis(4,3,j, 'SpacingHoriz', 0.03, 'SpacingVert',0.03,'Padding', 0.01, 'MarginRight', 0.1,...
         'MarginLeft',0.05,'MarginTop',0.01,'MarginBottom',0.05);
    imagesc(weightedPSM);
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
print('PSMS_dropout30_PY','-dpdf','-r300')



 

