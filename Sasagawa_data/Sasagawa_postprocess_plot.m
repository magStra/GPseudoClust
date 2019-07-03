%plot results without subsampling
%with subsampling
K = load('PSMs_mes_sub.mat');

clusterSolution = computeSummaryPSM_lmkk(K,2:12);

figure()
ax(4) = subaxis(1,4,4, 'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(clusterSolution.PSM(ind,ind));
set(gca,'YDir','normal')
title('lmkk');
colormap(hot);
ax(3) = subaxis(1,4,3, 'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
A = mean(K,3);
imagesc(A(ind,ind));
set(gca,'YDir','normal')
colormap(hot);
title('mean PSM');

%PY+PEAR
load('GPseudoClust/postProcPSM/mesPY.mat');
ax(1) = subaxis(1,4,1,'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(weightedPSM(ind,ind));
set(gca,'YDir','normal')
title('PY and PEAR');
ax(2) = subaxis(1,4,2,'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(weightedPSM_DP(ind,ind));
set(gca,'YDir','normal')
title('DPM and PEAR');
s2Pos = get(ax(4),'position');
s1Pos = get(ax(4),'position');
col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.01  s2Pos(2)  0.02  s1Pos(2)+s1Pos(3)*3.5])

set(gcf, 'PaperUnits', 'centimeters');
 x_width=17.6 ;y_width=4;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('PSMS_mes','-dpdf','-r350')

%second plot

A = zeros(600,600,4);
for j = 1:4
    A(:,:,j)=csvread(sprintf('PSM_mes%d.csv',j));
end
Amean = mean(A,3);
cl = sumClust(Amean,2:12);
[~,ind] = sort(cl);
load('GPseudoClust/postProcPSM/mesPY.mat');

set(0,'DefaultAxesFontSize', 7)
figure()
ax(1) = subaxis(2,2,1,'SpacingHoriz', 0.05, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.16,...
         'MarginLeft',0.07,'MarginTop',0.07,'MarginBottom',0.07);
imagesc(A(ind,ind,1));
set(gca,'YDir','normal')
colormap(hot);
title('all cells: chain 1','FontSize',7);
ax(2) = subaxis(2,2,2,'SpacingHoriz', 0.05, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.16,...
         'MarginLeft',0.07,'MarginTop',0.07,'MarginBottom',0.07);
imagesc(A(ind,ind,2));
set(gca,'YDir','normal')
colormap(hot);
title('all cells: chain 2','FontSize',7);
ax(4) = subaxis(2,2,4, 'SpacingHoriz', 0.05, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.16,...
         'MarginLeft',0.07,'MarginTop',0.07,'MarginBottom',0.07);
AA = mean(K,3);
imagesc(AA(ind,ind));
set(gca,'YDir','normal')
colormap(hot);
title('mean PSM','FontSize',7);
ax(3) = subaxis(2,2,3,'SpacingHoriz', 0.05, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.16,...
         'MarginLeft',0.07,'MarginTop',0.07,'MarginBottom',0.07);
imagesc(weightedPSM(ind,ind));
set(gca,'YDir','normal')
title('PY and PEAR','FontSize',7);  
s2Pos = get(ax(4),'position');
s1Pos = get(ax(4),'position');
col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.03  s2Pos(2)  0.05  s1Pos(2)+s1Pos(3)*2.25])
col.FontSize = 7;
set(gcf, 'PaperUnits', 'centimeters');
 x_width=8.8 ;y_width=8;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('PSMS_mesSummary','-dpdf','-r350')

%Identify groups with >80% pairwise co-clustering probability
load('GPseudoClust/postProcPSM/mesPY.mat');

mesGenes = importdata('mesGenes.csv');
mesGenes = mesGenes(2:end);
coSagasawa = coClustGeneLists(weightedPSM,mesGenes,0.8,0.01);

