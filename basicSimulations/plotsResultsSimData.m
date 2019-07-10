

%ARI first simulated data set
figure()
set(0,'DefaultAxesFontSize', 18)
set(0,'defaultfigurecolor',[1 1 1]);
ARI_sim1 = [ 0.94 ,0.8403228,1 ,1,1,1 , 0.59,0.9542601,1,1,1];
FMI_sim1 = [.95,0.8744034,1,1,1,1,0.74,0.9644013,1,1,1];
NMI_sim1 = [.95,0.9031884,1,1,1,1,0.78,0.9670273,1,1,1];
A = [ARI_sim1',FMI_sim1',NMI_sim1'];
A = importdata('clusterCompSim1.csv');
subplot(1,2,1);
barh(A.data');
xlim([0.2,1]);
legend('ARI','FMI','NMI','Location','NorthOutSide','Orientation','horizontal');
title('data set 1');
yticklabels({'mcl','SIMLR','PAM','hier.','SL+GCl','De+GCl','Mon 2','GPs+lmkk','GPs+PY','GPs+DP','GPs+mean'})
box off;

%second simulated data set
ARI_sim1 = [0.66,0.6962025,0.24,0.35,0.91,0.85,0.67,0.937847,1,1,1];
FMI_sim1 = [0.73,0.7589160,0.54,0.62,0.93,0.88,0.77,0.9516343,1,1,1];
NMI_sim1 = [0.81,0.7964319,0.42,0.58,0.9,0.88,0.75,0.94786641,1,1,1];
A = [ARI_sim1',FMI_sim1',NMI_sim1'];
A = importdata('clusterCompSim2.csv');
subplot(1,2,2);
barh(A.data');
xlim([0.2,1]);
legend('ARI','FMI','NMI','Location','NorthOutSide','Orientation','horizontal');
title('data set 2');
yticklabels({'mcl','SIMLR','PAM','hier.','SL+GCl','De+GCl','Mon 2','GPs+lmkk','GPs+PY','GPs+DP','GPs+mean'})
box off;
x_width=31 ;y_width=15;
set(gcf,'PaperUnits','centimeters');
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
set(gcf, 'PaperSize', [x_width y_width]);
print('simComp','-dpdf')


