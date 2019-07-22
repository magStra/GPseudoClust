function [] = simManyGenes(outputFile,inputSeed,nGenes)
rng(inputSeed);
set(0,'defaultFigureUnits','centimeters');
nCells = 60;
nClusters = randsample(2:5,1);
tau = (0.5:(nCells-0.5))/nCells;
[X, Y] = meshgrid(tau);
timeDiffs    = abs(X - Y);
%random cluster allocations
z = randsample(1:nClusters,nGenes,true);
sigmaW2 = rand(1,nClusters)*2+12;
sigmaE2 = randn(1,nClusters)+sigmaW2/3;
L = rand(1,nClusters)*0.4+0.3;
l = rand(1,nClusters)*0.2+0.9;
sigmaw2 = rand(1,nClusters)*0.2+0.4;
K = {};
Sigmas = {};
mu = {};
for j = 1:nClusters
    K{j} = sigmaW2(j)*(1+sqrt(3)*timeDiffs/L(j)).*exp(-sqrt(3)*timeDiffs/L(j));
    mu{j} = mvnrnd(zeros(1,60),K{j});
    Sigmas{j} = sigmaw2(j)*(1+sqrt(3)*timeDiffs/l(j)).*exp(-sqrt(3)*timeDiffs/l(j));
end

simData = zeros(nGenes,60);
for j = 1:nGenes
    simData(j,:)= mvnrnd(mu{z(j)},Sigmas{z(j)}+eye(nCells)*sigmaE2(z(j)));
end
csvwrite([outputFile,'.csv'],simData);
save([outputFile,'.mat'],'mu','nClusters','sigmaW2','sigmaE2','L','l','sigmaw2','z','nGenes');
end

