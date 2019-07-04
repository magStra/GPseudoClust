%combining the two cell lines
%run after running StumpfE.m, StumpfR.m, and StumpfPostprocessing.R

addpath(genpath('../'));
%add path to subaxis function
%https://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot
addpath(genpath('~/subaxis'));

K = zeros(94,94,192);
K1 = load('PSMs_StumpfE.mat');
K(:,:,1:96) = K1.PSMs;
K2 = load('PSMs_StumpfR.mat');
K(:,:,97:192) = K2.PSMs;


clusterSolution = computeSummaryPSM_lmkk(K,2:10);

figure()
ax(4) = subaxis(2,2,4,'SpacingHoriz', 0.05, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.13,...
         'MarginLeft',0.07,'MarginTop',0.07,'MarginBottom',0.07);
[x ind] = sort(clusterSolution.clustering);
imagesc(clusterSolution.PSM(ind,ind));
title('lmkk');
colormap(hot);
set(gca,'YDir','normal')

ax(3) = subaxis(2,2,3, 'SpacingHoriz', 0.05, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.13,...
         'MarginLeft',0.07,'MarginTop',0.07,'MarginBottom',0.07);
A = mean(K,3);
imagesc(A(ind,ind));
colormap(hot);
set(gca,'YDir','normal')

title('mean PSM');
%processing the results from PY
load('npSummaryPSMStumpf.mat');
ax(1) = subaxis(2,2,1,'SpacingHoriz', 0.05, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.13,...
         'MarginLeft',0.07,'MarginTop',0.07,'MarginBottom',0.07);
imagesc(npSummaryPSMStumpf.weightedPSM(ind,ind));
colormap(hot);
set(gca,'YDir','normal')

title('PY and PEAR');
ax(2) = subaxis(2,2,2,'SpacingHoriz', 0.05, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.13,...
         'MarginLeft',0.07,'MarginTop',0.07,'MarginBottom',0.07);
imagesc(npSummaryPSMStumpf.weightedPSM_DP(ind,ind));
colormap(hot);
set(gca,'YDir','normal')

title('DPM and PEAR');

s2Pos = get(ax(4),'position');
s1Pos = get(ax(4),'position');
col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.01  s2Pos(2)  0.05  s1Pos(2)+s1Pos(3)*2.15])

set(gcf, 'PaperUnits', 'centimeters');
 x_width=8.8 ;y_width=7;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);print('PSMS_Stumpf','-dpdf','-r350')

 %histogram of weights - PY
figure()
subplot(2,2,1);
histogram(npSummaryPSMStumpf.weights(1:96),'BinMethod','sqrt');
ylim([0,70]);
xlabel('weights')
title('PY: cell line 1');
subplot(2,2,2);
histogram(npSummaryPSMStumpf.weights(97:192),'BinMethod','sqrt');
ylim([0,70]);
xlabel('weights')
title('PY: cell line 2');
subplot(2,2,3);
histogram(npSummaryPSMStumpf.weights_DP(1:96),'BinMethod','sqrt');
ylim([0,70]);
xlabel('weights')
title('DP: cell line 1');
subplot(2,2,4);
histogram(npSummaryPSMStumpf.weights_DP(97:192),'BinMethod','sqrt');
ylim([0,70]);

xlabel('weights')
title('DP: cell line 2');

set(gcf, 'PaperUnits', 'centimeters');
 x_width=8.8 ;y_width=6;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);print('weights_Stumpf','-dpdf','-r350')

 
 %compare weights of two methods
 set(0,'DefaultAxesFontSize',6);
 figure()
 subplot(2,1,1);
 [xx, ind] = sort(npSummaryPSMStumpf.weights(1:96));
 scatter(1:96,xx,20,'MarkerFaceColor','g','MarkerEdgeColor','g','MarkerFaceAlpha',0.1);
 hold on;
 scatter(1:96,npSummaryPSMStumpf.weights_DP(ind),20,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',0.1);
 %leg=legend('PY+PEAR','DPM+PEAR','Location','NorthOutside');
 xlabel('ascending order for PY and PEAR','FontSize',8);
 ylabel('weight','FontSize',8);
 title('Cell line 1','FontSize',8);
 %leg.FontSize = 6;

 subplot(2,1,2);
 [xx, ind] = sort(npSummaryPSMStumpf.weights(97:192));
 scatter(1:96,xx,20,'MarkerFaceColor','g','MarkerEdgeColor','g','MarkerFaceAlpha',0.1);
 hold on;
 scatter(1:96,npSummaryPSMStumpf.weights_DP(ind+96),20,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',0.1);
 %leg=legend('PY+PEAR','DPM+PEAR','Location','NorthOutside');
 %leg.FontSize = 6;
  xlabel('ascending order for PY and PEAR','FontSize',8);
 ylabel('weight','FontSize',8);

  title('Cell line 2','FontSize',8);

% leg.FontSize = 6;
set(gcf, 'PaperUnits', 'centimeters');
 x_width=8.8 ;y_width=9;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);print('weights_Stumpf1','-dpdf','-r350')

 
[co,noCo] = coClustGeneLists(npSummaryPSMStumpf.weightedPSM,importdata('StumpfGeneNames.csv'),0.8,0.01);

 