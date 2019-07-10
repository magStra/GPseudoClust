close all;
clear all;
addpath(genpath('~/subaxis'));

%plotting first simulated data set
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)
set(0,'defaultfigurecolor',[1 1 1]);

figure()
%subaxis(1,2,1, 'SpacingHoriz', 0.05, 'SpacingVert',0.15,'Padding', 0.01, 'MarginRight', 0.02,...
%         'MarginLeft',0.1,'MarginTop',0.1,'MarginBottom',0.2);
subplot(1,2,1);
simData = csvread('simDataClust.csv');
tau = (1:60)/60;
c1 = plot(tau,simData(1,:),'go','MarkerEdgeColor','g','MarkerSize',1);
hold on;
plot(tau,simData(1:8,:),'go','MarkerEdgeColor','g','MarkerSize',1);
hold on;
c2 = plot(tau,simData(9,:),'mo','MarkerEdgeColor','m','MarkerSize',1);
hold on;
plot(tau,simData(9:12,:),'mo','MarkerEdgeColor','m','MarkerSize',1);
hold on;
c3 =plot(tau,simData(13,:),'ko','MarkerEdgeColor','k','MarkerSize',1);
hold on;
plot(tau,simData(13:24,:),'ko','MarkerEdgeColor','k','MarkerSize',1);
hold on;
c4 = plot(tau,simData(25,:),'bo','MarkerEdgeColor','b','MarkerSize',1);
hold on;
plot(tau,simData(25:40,:),'bo','MarkerEdgeColor','b','MarkerSize',1);
hold on;
c5 = plot(tau,simData(41,:),'co','MarkerEdgeColor','c','MarkerSize',1);
plot(tau,simData(41:52,:),'co','MarkerEdgeColor','c','MarkerSize',1);
xlabel('pseudotime','FontSize',12);
ylabel('expression level','FontSize',12);
title('Data set 1','FontSize',12);
set(gca,'box','off');
ax=gca;
ylim([-10,10]);

% leg = legend([c1,c2,c3,c4,c5],'cluster 1','cluster 2','cluster 3','cluster 4','cluster 5');
% leg.FontSize = 24;

%subaxis(1,2,2, 'SpacingHoriz', 0.1, 'SpacingVert',0.15,'Padding', 0.01, 'MarginRight', 0.02,...
 %        'MarginLeft',0.1,'MarginTop',0.1,'MarginBottom',0.2);
 subplot(1,2,2);
simData = csvread('simDataClust2.csv');
tau = (1:60)/60;
c1 = plot(tau,simData(1,:),'go','MarkerEdgeColor','g','MarkerSize',1);
hold on;
plot(tau,simData(1:8,:),'go','MarkerEdgeColor','g','MarkerSize',1);
hold on;
c2 = plot(tau,simData(9,:),'mo','MarkerEdgeColor','m','MarkerSize',1);
hold on;
plot(tau,simData(9:12,:),'mo','MarkerEdgeColor','m','MarkerSize',1);
hold on;
c3 =plot(tau,simData(13,:),'ko','MarkerEdgeColor','k','MarkerSize',1);
hold on;
plot(tau,simData(13:24,:),'ko','MarkerEdgeColor','k','MarkerSize',1);
hold on;
c4 = plot(tau,simData(25,:),'bo','MarkerEdgeColor','b','MarkerSize',1);
hold on;
plot(tau,simData(25:40,:),'bo','MarkerEdgeColor','b','MarkerSize',1);
hold on;
c5 = plot(tau,simData(41,:),'co','MarkerEdgeColor','c','MarkerSize',1);
plot(tau,simData(41:52,:),'co','MarkerEdgeColor','c','MarkerSize',1);
xlabel('pseudotime','FontSize',12);
ylabel('expression level','FontSize',12);
title('Data set 2','FontSize',12);
set(gca,'box','off');
ax=gca;
ylim([-10,10]);
set(gcf, 'PaperUnits', 'centimeters');
x_width=17.6 ;y_width=6;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
set(gcf, 'PaperSize', [x_width y_width]);
print('simData','-dpdf');

set(0,'DefaultAxesFontSize', 10)

