addpath(genpath('../'));
%add path to subaxis function
%https://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot
addpath(genpath('~/subaxis'));

%compute summary matrix representation using lmkk on all 96 chains
K = load('PSMsShalek.mat');
K = K.PSMsShalek;
clusterSolution = computeSummaryPSM_lmkk(K,2:12);

set(0,'defaultAxesFontSize',6);


figure()
load('npSummaryPSMShalek.mat');

%ordering of the rows and columns of the matrices to make the cluster
%structure visible. Please note that the labels in the summary clustering
%may swap across runs, and therefore the blocks of genes with high
%co-clustering probabilities may be permuted, and the plot therefore look
%slightly different for each run (blocks differently arranged), even if the PSMs themselves are almost
%identical. 
[x ind] = sort(npSummaryPSMShalek.sumClustPEAR);

ax(1) = subaxis(1,4,1,'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(npSummaryPSMShalek.weightedPSM(ind,ind));
set(gca,'yDir','normal');
title('PY and PEAR');

ax(2) = subaxis(1,4,2,'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(npSummaryPSMShalek.weightedPSM_DP(ind,ind));
set(gca,'yDir','normal');
title('DPM and PEAR');

A = mean(K,3);
ax(3) = subaxis(1,4,3,'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(A(ind,ind));
set(gca,'yDir','normal');
title('mean PSM');

ax(4) = subaxis(1,4,4, 'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(clusterSolution.PSM(ind,ind));
set(gca,'yDir','normal');
title('lmkk');

s2Pos = get(ax(4),'position');
s1Pos = get(ax(4),'position');
col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.01  s2Pos(2)  0.02  s1Pos(2)+s1Pos(3)*3.45])

set(gcf, 'PaperUnits', 'centimeters');
 x_width=17.6 ;y_width=4.5;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('PSMS_Shalek','-dpdf','-r350')


%checking stability with regard to number of subsampled chains

%24 subsampled chains
K = load('PSMsShalek.mat');
K = K.PSMsShalek(:,:,1:24);
clusterSolution24 = computeSummaryPSM_lmkk(K,2:12);
figure()
load('npSummaryPSMShalek24.mat');
[x ind] = sort(npSummaryPSMShalek.sumClustPEAR);
ax(1) = subaxis(1,4,1,'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(npSummaryPSMShalek24.weightedPSM(ind,ind));
set(gca,'yDir','normal');
title('PY and PEAR');

ax(2) = subaxis(1,4,2,'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(npSummaryPSMShalek24.weightedPSM_DP(ind,ind));
set(gca,'yDir','normal');
title('DPM and PEAR');

A = mean(K,3);
ax(3) = subaxis(1,4,3,'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(A(ind,ind));
set(gca,'yDir','normal');
title('mean PSM');

ax(4) = subaxis(1,4,4, 'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(clusterSolution24.PSM(ind,ind));
set(gca,'yDir','normal');
title('lmkk');

s2Pos = get(ax(4),'position');
s1Pos = get(ax(4),'position');
col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.01  s2Pos(2)  0.02  s1Pos(2)+s1Pos(3)*3.45])

set(gcf, 'PaperUnits', 'centimeters');
 x_width=17.6 ;y_width=4.5;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('PSMS_Shalek24','-dpdf','-r350')

%24 subsampled chains
K = load('PSMsShalek.mat');
K = K.PSMsShalek(:,:,1:24);
clusterSolution4 = computeSummaryPSM_lmkk(K,2:12);

figure()
load('npSummaryPSMShalek4.mat');
[x ind] = sort(npSummaryPSMShalek.sumClustPEAR);
ax(1) = subaxis(1,4,1,'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(npSummaryPSMShalek4.weightedPSM(ind,ind));
set(gca,'yDir','normal');
title('PY and PEAR');

ax(2) = subaxis(1,4,2,'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(npSummaryPSMShalek4.weightedPSM_DP(ind,ind));
set(gca,'yDir','normal');
title('DPM and PEAR');

A = mean(K,3);
ax(3) = subaxis(1,4,3,'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(A(ind,ind));
set(gca,'yDir','normal');
title('mean PSM');

ax(4) = subaxis(1,4,4, 'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(clusterSolution4.PSM(ind,ind));
set(gca,'yDir','normal');
title('lmkk');

s2Pos = get(ax(4),'position');
s1Pos = get(ax(4),'position');
col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.01  s2Pos(2)  0.02  s1Pos(2)+s1Pos(3)*3.45])

set(gcf, 'PaperUnits', 'centimeters');
 x_width=17.6 ;y_width=4.5;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('PSMS_Shalek4','-dpdf','-r350')



%for PY and 96 chains, compute lists of genes co-clustered with high
%posterior probability (80%)
load('npSummaryPSMShalek.mat');
genesShalek = importdata('genesShalek.csv');
genesShalek = [genesShalek(2:end)];
coShalek = coClustGeneLists(npSummaryPSMShalek.weightedPSM,genesShalek,0.8,0.01);
save('Shalek_coClustering08','coShalek');

%assignements of genes to the 80% groups
indices = {};
for j = 1:length(coShalek)
    indices{j} = [];
    for k = 1:length(coShalek{j})
        [~, idx] = ismember(coShalek{j}{k}, genesShalek);
        indices{j} = [indices{j},idx];
    end
end
