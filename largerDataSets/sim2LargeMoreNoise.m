function [] = sim2LargeMoreNoise(outputFile,inputSeed)
rng(inputSeed);
set(0,'defaultFigureUnits','centimeters');
nGenes = randsample(20:30,1);
nCells = 9000;
nClusters = randsample(2:5,1);
tau = (0.5:(nCells-0.5))/nCells;
[X, Y] = meshgrid(tau);
timeDiffs    = abs(X - Y);
%random cluster allocations
z = randsample(1:nClusters,nGenes,true);
sigmaW2 = rand(1,nClusters)*2+12;
sigmaE2 = randn(1,nClusters)+sigmaW2;
L = rand(1,nClusters)*0.4+0.3;
l = rand(1,nClusters)*0.2+0.9;
sigmaw2 = rand(1,nClusters)*0.2+0.4;
K = {};
Sigmas = {};
mu = {};
for j = 1:nClusters
    K{j} = sigmaW2(j)*(1+sqrt(3)*timeDiffs/L(j)).*exp(-sqrt(3)*timeDiffs/L(j));
    mu{j} = mvnrnd(zeros(1,9000),K{j});
    Sigmas{j} = sigmaw2(j)*(1+sqrt(3)*timeDiffs/l(j)).*exp(-sqrt(3)*timeDiffs/l(j));
end

simData = zeros(nGenes,9000);
for j = 1:nGenes
    simData(j,:)= mvnrnd(mu{z(j)},Sigmas{z(j)}+eye(nCells)*sigmaE2(z(j)));
end
colors = ['g','m','b','r','k'];
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultaxesfontsize',12)
figure()
for j = 1:nGenes
    plot(tau,simData(j,:),'o','Color',colors(z(j)),'MarkerSize',3);
    hold on;
end
 set(gcf, 'PaperPosition', [0 0 8.8 6]);
 set(gcf,'PaperSize',[8.8 6]);
 print([outputFile,'.pdf'],'-dpdf');
csvwrite([outputFile,'.csv'],simData);
save([outputFile,'.mat'],'mu','nClusters','sigmaW2','sigmaE2','L','l','sigmaw2','z','nGenes');
close all;
figure()
for j = 1:nGenes
    plot(tau,simData(j,:),'o','Color',colors(z(j)),'MarkerSize',3);
    hold on;
end
set(gcf, 'PaperPosition', [0 0 8.8 6]);
set(gcf,'PaperSize',[8.8 6]);
print([outputFile,'.pdf'],'-dpdf');
csvwrite([outputFile,'.csv'],simData);
close all;
end
%code to read in samples of concentration parameters from output files
% alphas = [];
% for k = 1:5
%     for j = 1:48
%         File = sprintf('simLargeMoreNoise%d_Results_Chain%d.csv',k,j);
%         xx = dlmread(File,',',[0 0 2500-1 0]);
%         alphas = [alphas xx];
%     end
%     csvwrite(sprintf('alphasSimLargeMoreNoise%d.csv',k),alphas);
% end