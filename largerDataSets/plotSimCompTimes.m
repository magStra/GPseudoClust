figure()
load('simLargeTime1.mat');
boxplot(elapsedTime/60,'Labels',{'10/30','20/60','30/90','40/120','50/150'});
xlabel('cells per capture time/total number of cells');
ylabel('elapsed time (min) for 5000 thinned samples'); 
set(gcf, 'PaperUnits', 'centimeters');
 x_width=8.8 ;y_width=5;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('compTimes5000','-dpdf','-r350')
