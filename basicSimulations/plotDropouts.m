%postprocessing

set(0,'DefaultAxesFontSize', 12)
set(0,'defaultfigurecolor',[1 1 1]);
addpath(genpath(pwd));
close all

figure()

A= csvread('scores_dropout_lmkk.csv')';
subplot(2,1,1);
boxplot(A,'Labels',{'ARI','FMI','NMI'});
ylim([0.4,1]);
set(gcf, 'PaperUnits', 'centimeters');
title('lower dropout levels','FontSize',12);

A= csvread('scores_dropout30_lmkk.csv')';
subplot(2,1,2);
boxplot(A,'Labels',{'ARI','FMI','NMI'});
ylim([0.4,1]);
set(gcf, 'PaperUnits', 'centimeters');
title('higher dropout levels','FontSize',12);

x_width=8.8 ;y_width=12;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('scores_Dropout_lmkk','-dpdf','-r350')


%the same with consensus clustering

figure()

A= csvread('scores_dropout_consensus.csv')';
subplot(2,1,1);
boxplot(A,'Labels',{'ARI','FMI','NMI'});
ylim([0.4,1]);
set(gcf, 'PaperUnits', 'centimeters');
title('lower dropout levels','FontSize',12);

A= csvread('scores_dropout30_consensus.csv')';
subplot(2,1,2);
boxplot(A,'Labels',{'ARI','FMI','NMI'});
ylim([0.4,1]);set(gcf, 'PaperUnits', 'centimeters');
title('higher dropout levels','FontSize',12);

x_width=8.8 ;y_width=12;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('scores_Dropout_consensus','-dpdf','-r350')

%same with Pitman-Yor

figure()

A= csvread('scores_dropout_PitmanYor.csv')';
subplot(2,1,1);
boxplot(A,'Labels',{'ARI','FMI','NMI'});
ylim([0.4,1]);
set(gcf, 'PaperUnits', 'centimeters');
title('lower dropout levels','FontSize',12);

A= csvread('scores_dropout30_PitmanYor.csv')';
subplot(2,1,2);
boxplot(A,'Labels',{'ARI','FMI','NMI'});
ylim([0.4,1]);
set(gcf, 'PaperUnits', 'centimeters');
title('higher dropout levels','FontSize',12);

x_width=8.8 ;y_width=12;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('scores_Dropout_PitManYor','-dpdf','-r350')

figure()

A= csvread('scores_dropout_DPM.csv')';
subplot(2,1,1);
boxplot(A,'Labels',{'ARI','FMI','NMI'});
ylim([0.4,1]);
set(gcf, 'PaperUnits', 'centimeters');
title('lower dropout levels','FontSize',12);

A= csvread('scores_dropout30_DPM.csv')';
subplot(2,1,2);
boxplot(A,'Labels',{'ARI','FMI','NMI'});
ylim([0.4,1]);
set(gcf, 'PaperUnits', 'centimeters');
title('higher dropout levels','FontSize',12);

x_width=8.8 ;y_width=12;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('scores_Dropout_DPM','-dpdf','-r350')




figure()

A= csvread('scores_dropoutMixed_consensus.csv')';
subplot(2,2,1);
boxplot(A,'Labels',{'ARI','FMI','NMI'});
ylim([0.4,1]);
title('mean PSM','FontSize',12);
subplot(2,2,2);
A= csvread('scores_dropoutMixed_lmkk.csv')';
boxplot(A,'Labels',{'ARI','FMI','NMI'});
ylim([0.4,1]);title('lmkk');
A= csvread('scores_dropoutMixed_PitmanYor.csv')';
subplot(2,2,3);
boxplot(A,'Labels',{'ARI','FMI','NMI'});
ylim([0.4,1]);
title('PY + PEAR','FontSize',12);
subplot(2,2,4);
A= csvread('scores_dropoutMixed_DPM.csv')';
boxplot(A,'Labels',{'ARI','FMI','NMI'});
ylim([0.4,1]);title('DPM + PEAR','FontSize',12);
x_width=8.8 ;y_width=6;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('scores_Dropout_mixed','-dpdf','-r350')


set(0,'DefaultAxesFontSize', 23)

%different figure of same thing
figure()

A= csvread('scores_dropoutMixed_consensus.csv')';
B= csvread('scores_dropoutMixed_lmkk.csv')';
C= csvread('scores_dropoutMixed_PitmanYor.csv')';
D= csvread('scores_dropoutMixed_DPM.csv')';

subaxis(1,2,1,'SpacingHoriz', 0.05, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.03,...
         'MarginLeft',0.03,'MarginTop',0.1,'MarginBottom',0.1);
boxplot([C(:,1),D(:,1),A(:,1),B(:,1)],'Labels',{'PY+PEAR','DPM+PEAR','mean psm','lmkk'});
ylim([0.4,1]);
title('ARI','FontSize',24);
subaxis(1,2,2,'SpacingHoriz', 0.05, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.03,...
         'MarginLeft',0.03,'MarginTop',0.1,'MarginBottom',0.1);
boxplot([C(:,3),D(:,3),A(:,3),B(:,3)],'Labels',{'PY+PEAR','DPM+PEAR','mean psm','lmkk'});
ylim([0.4,1]);
title('NMI','FontSize',24);
x_width=17.6 ;y_width=5;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('scores_Dropout_mixed1','-dpdf','-r350')



