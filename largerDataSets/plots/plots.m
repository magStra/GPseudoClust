addpath(genpath('../'));
addpath(genpath('~/Documents/GPseudoRankCode'));%add path with subaxis function

nCells = 9000;
tau = (0.5:(nCells-0.5))/nCells;
colors = ['g','m','b','r','k'];
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultaxesfontsize',16);
set(0, 'defaultFigureUnits', 'centimeters');
for j = 1:10
    j
    figure()
    for k = 1:3
        ax(k+1) = subaxis(1,4,k+1,'SpacingHoriz', 0.02, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.03,'MarginTop',0.07,'MarginBottom',0.07);
        colormap(hot);
        A = csvread(sprintf('../PSMs_PY_%d_%d.csv',10*k,j));
        imagesc(A);
        set(gca,'YDir','normal');
        title(sprintf('%d subsampled cells',k*10),'FontSize',16);
    end
    simData = csvread(sprintf('../simLarge%d.csv',j));
    load(sprintf('../simLarge%d.mat',j));
    ax(1) = subaxis(1,4,1,'SpacingHoriz', 0.02, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.03,'MarginTop',0.07,'MarginBottom',0.07);
    for jj = 1:nGenes
        plot(tau,simData(jj,:),'o','Color',colors(z(jj)),'MarkerSize',0.3);
        hold on;
    end
    title('data set','FontSize',16);
    s2Pos = get(ax(4),'position');
    s1Pos = get(ax(4),'position');
    col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.01  s2Pos(2)  0.02  s1Pos(2)+s1Pos(3)*4.05]);
    set(gcf, 'PaperPosition', [0 0 17.6,4]);
    set(gcf,'PaperSize',[17.6,4]);

    print(sprintf('simLarge%d',j),'-dpdf','-r350')

end
close all;

xx = [2,4,8,12,48,96];
for j = 1:10
    figure()
    for k = 1:6
        ax(k+1) = subaxis(2,4,k+1,'SpacingHoriz', 0.02, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.03,'MarginTop',0.07,'MarginBottom',0.07);
        colormap(hot);
        A = csvread(sprintf('../PSMs_PY%d_%d.csv',xx(k),j));
        imagesc(A);
        set(gca,'YDir','normal');
        title(sprintf('%d chains',xx(k)),'FontSize',16);
    end
    simData = csvread(sprintf('../simLarge%d.csv',j));
    load(sprintf('../simLarge%d.mat',j));
    ax(1) = subaxis(2,4,1,'SpacingHoriz', 0.02, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.03,'MarginTop',0.07,'MarginBottom',0.07);
    for jj = 1:nGenes
        plot(tau,simData(jj,:),'o','Color',colors(z(jj)),'MarkerSize',0.3);
        hold on;
    end
    title('data set','FontSize',16);
    s2Pos = get(ax(7),'position');
    s1Pos = get(ax(7),'position');
    col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.04  s2Pos(2)  0.04  s1Pos(2)+s1Pos(3)*1.45]);
    set(gcf, 'PaperPosition', [0 0 17.6,8]);
    set(gcf,'PaperSize',[17.6,8]);

    print(sprintf('simLarge10_%d',j),'-dpdf','-r350')

end
close all;

for j = 1:5
    figure()
    for k = 1:3
        ax(k+1) = subaxis(1,4,k+1,'SpacingHoriz', 0.02, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.03,'MarginTop',0.07,'MarginBottom',0.07);
        colormap(hot);
        A = csvread(sprintf('../PSMs_PY_%d_MoreNoise%d.csv',10*k,j));
        imagesc(A);
        set(gca,'YDir','normal');
        title(sprintf('%d subsampled cells',k*10),'FontSize',16);
    end
    simData = csvread(sprintf('../simLargeMoreNoise%d.csv',j));
    load(sprintf('../simLargeMoreNoise%d.mat',j));
    ax(1) = subaxis(1,4,1,'SpacingHoriz', 0.02, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.03,'MarginTop',0.07,'MarginBottom',0.07);
    for jj = 1:nGenes
        plot(tau,simData(jj,:),'o','Color',colors(z(jj)),'MarkerSize',0.3);
        hold on;
    end
    title('data set','FontSize',16);
    s2Pos = get(ax(4),'position');
    s1Pos = get(ax(4),'position');
    col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.01  s2Pos(2)  0.02  s1Pos(2)+s1Pos(3)*4.05]);
    set(gcf, 'PaperPosition', [0 0 17.6,4]);
    set(gcf,'PaperSize',[17.6,4]);

    print(sprintf('simLargeMoreNoise%d',j),'-dpdf','-r350')

end
close all;

xx = [2,4,8,12,48,96];
for j = 1:5
    figure()
    for k = 1:6
        ax(k+1) = subaxis(2,4,k+1,'SpacingHoriz', 0.02, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.03,'MarginTop',0.07,'MarginBottom',0.07);
        colormap(hot);
        A = csvread(sprintf('../PSMs_PY%d_MoreNoise%d.csv',xx(k),j));
        imagesc(A);
        set(gca,'YDir','normal');
        title(sprintf('%d chains',xx(k)),'FontSize',16);
    end
    simData = csvread(sprintf('../simLargeMoreNoise%d.csv',j));
    load(sprintf('../simLargeMoreNoise%d.mat',j));
    ax(1) = subaxis(2,4,1,'SpacingHoriz', 0.02, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.03,'MarginTop',0.07,'MarginBottom',0.07);
    for jj = 1:nGenes
        plot(tau,simData(jj,:),'o','Color',colors(z(jj)),'MarkerSize',0.3);
        hold on;
    end
    title('data set','FontSize',16);
    s2Pos = get(ax(7),'position');
    s1Pos = get(ax(7),'position');
    col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.04  s2Pos(2)  0.04  s1Pos(2)+s1Pos(3)*1.45]);
    set(gcf, 'PaperPosition', [0 0 17.6,8]);
    set(gcf,'PaperSize',[17.6,8]);

    print(sprintf('simLarge10_MoreNoise%d',j),'-dpdf','-r350')

end
close all;