figure();
A=csvread('Shalek_ARI.csv');
imagesc(A);
colormap(gray)
colorbar;
col = colorbar;
col.FontSize = 14;
set(0,'defaultAxesFontSize',16)
x_width=8.8 ;y_width=4;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
set(gcf, 'PaperSize', [x_width y_width]);
xticklabels({'mcl','SIMLR','PAM','hier.','SL+GCl','De+GCl','Mon 2'})
yticklabels({'mcl','SIMLR','PAM','hier.','SL+GCl','De+GCl','Mon 2'})
print('Rand_Shalek','-dpdf')

figure()
A=csvread('Shalek_ARIExt.csv');
imagesc(A);
colorbar;
col = colorbar;
col.FontSize = 10;
set(0,'defaultAxesFontSize',13)
x_width=8.8 ;y_width=4;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
set(gcf, 'PaperSize', [x_width y_width]);
xticks(1:size(A,1));
xticklabels({'mcl','SIMLR','PAM','hier.','SL+GCl','De+GCl','Mon 2','sample 1','sample 2','sample 3','sample 4',...
    'sample 5','sample 6','sample 7','sample 8'})
xtickangle(45)
yticks(1:size(A,1));
yticklabels({'mcl','SIMLR','PAM','hier.','SL+GCl','De+GCl','Mon 2','sample 1','sample 2','sample 3','sample 4',...
    'sample 5','sample 6','sample 7','sample 8'})
set(gca,'yDir','normal');
print('Rand_ShalekExt','-dpdf')
