% clear all
%To identify the branches using diffusion pseudotime (Haghverdi et al.
%2016), download the code provided as supplementary materials to 
%Haghverdi et al. 2016
%https://www.nature.com/articles/nmeth.3971#supplementary-information

pathToDiffusionPseudotime  = '~/Documents/GPseudoRankCode/DPT_inMatlab/dpt/examples/';%edit!
addpath(genpath(pathToDiffusionPseudotime));
run ESC_qpcr 
for j = 1:3
   data1 = data(Branch==j,:);
   CT1 = labels(Branch==j);
   csvwrite(sprintf('ESC_%d.csv',j),data1');
   csvwrite(sprintf('ESC_CT_%d.csv',j),CT1);
end

%now run GPseudoClust on each of the branches
%(ESC1.m,ESC2.m,ESC3.m) and compute the summary PSMs (npSummaryPSM_ESC.R)

addpath(genpath('../'));

load('npSummaryPSM_ESC1.mat');
weightedPSM1 = npSummaryPSM1.weightedPSM;
load('npSummaryPSM_ESC2.mat');
weightedPSM2 = npSummaryPSM2.weightedPSM;
load('npSummaryPSM_ESC3.mat');
weightedPSM3 = npSummaryPSM3.weightedPSM;

set(0,'DefaultAxesFontSize', 7)
%add path to subaxis function
%https://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot
addpath(genpath('~/subaxis'));


figure()
ax(1) = subaxis(2,2,1,'SpacingHoriz', 0.05, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.13,...
         'MarginLeft',0.07,'MarginTop',0.07,'MarginBottom',0.07);

meanPSM = 1/3*(weightedPSM1+weightedPSM2+weightedPSM3);
sumCl = sumClust(meanPSM,4:10);
[x ind] = sort(sumCl);
imagesc(weightedPSM1(ind,ind));
set(gca,'YDir','normal')

colormap(gray);
title('trunk','FontSize',8);
ax(2) = subaxis(2,2,3,'SpacingHoriz', 0.05, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.13,...
         'MarginLeft',0.07,'MarginTop',0.07,'MarginBottom',0.07);

imagesc(weightedPSM2(ind,ind));
set(gca,'YDir','normal')

colormap(gray);
title('erythroid branch','FontSize',8);
ax(3) = subaxis(2,2,2,'SpacingHoriz', 0.05, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.13,...
         'MarginLeft',0.07,'MarginTop',0.07,'MarginBottom',0.07);
imagesc(weightedPSM3(ind,ind));
set(gca,'YDir','normal')

colormap(gray);
title('endothelial branch','FontSize',8);
ax(4) = subaxis(2,2,4,'SpacingHoriz', 0.05, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.13,...
         'MarginLeft',0.07,'MarginTop',0.07,'MarginBottom',0.07);

imagesc(meanPSM(ind,ind));
set(gca,'YDir','normal')

colormap(gray);
title('mean PSM','FontSize',8);
s2Pos = get(ax(4),'position');
s1Pos = get(ax(4),'position');
col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.01  s2Pos(2)  0.05  s1Pos(2)+s1Pos(3)*2.17])

set(gcf, 'PaperUnits', 'centimeters');
 x_width=8.8 ;y_width=7;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('PSMS_ESC_PY','-dpdf','-r350')

%for each of the PSMs compute adjacency matrix, where there is an edge if
%the probability of co-clustering is above 80%
weightedPSM1a =  weightedPSM1 > 0.8;
figure()
plot(graph(weightedPSM1a)); 
box off;
axis off;
 x_width=4 ;y_width=2.5;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('graph_ESC_root','-dpdf','-r350')


weightedPSM2a =  weightedPSM2 > 0.8;
figure()
plot(graph(weightedPSM2a)); 
box off;
axis off;
 x_width=4 ;y_width=2.5;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('graph_ESC_branch2','-dpdf','-r350')

weightedPSM3a =  weightedPSM3 > 0.8;
figure()
plot(graph(weightedPSM3a)); 
box off;
axis off;
 x_width=4 ;y_width=2.5;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('graph_ESC_branch3','-dpdf','-r350')

meanPSMa =  meanPSM > 0.8;
figure()
plot(graph(meanPSMa)); 
box off;
axis off;
 x_width=4 ;y_width=2.5;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('graph_ESC_mean','-dpdf','-r350')

%groups of genes with posterior pairwise co-clustering probability > 0.8
load ESC_qpcr_Goettgens.mat;
[co,noCo] = coClustGeneLists(weightedPSM1,Genes_analysed,0.8,0.01);

[co2,noCo2] = coClustGeneLists(weightedPSM2,Genes_analysed,0.8,0.01);

[co3,noCo3] = coClustGeneLists(weightedPSM3,Genes_analysed,0.8,0.01);



%plotting PSMs for each branch separately
load ESC_qpcr_Goettgens.mat;

%root
[xx ind] = sort(npSummaryPSM1.sumClustPEAR);
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultaxesfontsize',20)
fig = figure()
imagesc(weightedPSM1(ind,ind))
colormap(gray);
set(gca,'YDir','normal')
xticklabels(Genes_analysed(ind));
xticks(1:42);
xtickangle(90);
yticklabels(Genes_analysed(ind));
yticks(1:42);
x_width=17 ;y_width=12;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('PSM_ESC_trunk','-dpng')


fig = figure()
[xx ind] = sort(npSummaryPSM2.sumClustPEAR);
imagesc(weightedPSM2(ind,ind))
colormap(gray);
set(gca,'YDir','normal')
xticklabels(Genes_analysed(ind));
xticks(1:42);
xtickangle(90);
yticklabels(Genes_analysed(ind));
yticks(1:42);
x_width=17 ;y_width=12;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('Erythroid branch','-dpng')


fig = figure()
[xx ind] = sort(npSummaryPSM3.sumClustPEAR);
imagesc(weightedPSM3(ind,ind))
colormap(gray);
set(gca,'YDir','normal')
xticklabels(Genes_analysed(ind));
xticks(1:42);
xtickangle(90);
yticklabels(Genes_analysed(ind));
yticks(1:42);
x_width=17 ;y_width=12;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('Endothelial branch','-dpng')

