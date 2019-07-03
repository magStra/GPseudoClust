function [mu,nClusters,b,c,sigma2,z,nGenes] = generateDataLinear(fileName)
%simulates data from clusters of GPs with linear covariance matrices
%output: 
set(0,'defaultFigureUnits','centimeters');
nGenes = randsample(20:30,1);
nCells = 60;
nClusters = randsample(2:5,1);
tau = (0.5:(nCells-0.5))/nCells;
[X, Y] = meshgrid(tau);
%random cluster allocations
z = randsample(1:nClusters,nGenes,true);
b = zeros(1,nClusters);
c = zeros(1,nClusters);
sigma2 = zeros(1,nClusters);
for j = 1:nClusters
    b(j) = (nClusters - 2*j)*3 + randn();
    c(j) = randn()+j*2;
    sigma2(j) = 1+rand()*0.2;
end
mu = zeros(nClusters,nCells);
for j = 1:nClusters
    mu(j,:) = b(j).*tau+c(j);
end
 simData = zeros(nGenes,nCells);
for j = 1:nGenes
    simData(j,:)= mvnrnd(mu(z(j),:),eye(nCells)*sigma2(z(j)));
end
colors = ['g','m','b','r','k'];
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultaxesfontsize',12)
for j = 1:nGenes
   plot(tau,simData(j,:),'o','Color',colors(z(j)),'MarkerSize',3);
   hold on;
end
 set(gcf, 'PaperPosition', [0 0 8.8 6]);
 set(gcf,'PaperSize',[8.8 6]);
 print([fileName,'.pdf'],'-dpdf');
csvwrite([fileName,'.csv'],simData);
save([fileName,'.mat'],'mu','nClusters','c','b','sigma2','z','nGenes');
close all;
end

