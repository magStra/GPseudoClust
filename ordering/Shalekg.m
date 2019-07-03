%Using GPseudoClustPlus to obtain samples from the joint distribution of 
%orders and cluster allocations. 

set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 14)
set(0,'defaultfigurecolor',[1 1 1]);
addpath(genpath('../'));

%Pre-order genes 

A = csvread('Shalek13.csv');
B = csvread('Shalek13AV.csv');%mean expression of antiviral genes as in Shalek
% et al, 2013. 
[xx, ind] = sort(B,'ascend');
csvwrite('Shalek13PreOrderedAV.csv',A(:,ind));
fHandle          =  @GPseudoClustPlus;  
fileName         = 'Shalek13PreOrderedAV.csv';
A = csvread('Shalek13PreOrderedAV.csv');
[nGenes,nCells]  = size(A);
nSamples         = 20000;  
verbose          = false; 
inputSeed        = NaN;    
permuteData      = true;
b = 0.01;
captureTimes     = [repmat(1,1,3),repmat(2,1,3),repmat(3,1,3),...
    repmat(4,1,3),repmat(5,1,3),repmat(6,1,3)];
adjustForCellSize = false;
pp               = [ 0 0 0.5 0.5];
delta            = [1/8000,1/8000];
n0 = 4;
n3 = 3;
n3a = 5;
kk = 0.01;
jj = 0.01;

parpool(12);
parfor j = 301:336
 tic
 feval(fHandle, fileName, j, nSamples, verbose,inputSeed,...
       permuteData,repmat(1,1,nCells),b,adjustForCellSize,1:nCells,...
       n0,n3,n3a,jj,kk,pp,delta);
 toc% 
 end
%
% Elapsed time is 4626.300634 seconds.
% Elapsed time is 4888.525983 seconds.
% Elapsed time is 5001.064181 seconds.
% Elapsed time is 5430.624619 seconds.
% Elapsed time is 5522.823121 seconds.
% Elapsed time is 5629.031060 seconds.
% Elapsed time is 5686.836983 seconds.
% Elapsed time is 5697.228781 seconds.
% Elapsed time is 5775.514154 seconds.
% Elapsed time is 5777.833345 seconds.
% Elapsed time is 5818.694876 seconds.
% Elapsed time is 5856.220134 seconds.
% Elapsed time is 5487.011220 seconds.
% Elapsed time is 5251.232889 seconds.
% Elapsed time is 5497.733604 seconds.
% Elapsed time is 5694.623610 seconds.
% Elapsed time is 5649.478740 seconds.
% Elapsed time is 5566.005927 seconds.
% Elapsed time is 5571.943266 seconds.
% Elapsed time is 5705.395601 seconds.
% Elapsed time is 5776.357749 seconds.
% Elapsed time is 5729.238302 seconds.
% Elapsed time is 5690.881285 seconds.
% Elapsed time is 5732.525444 seconds.
% Elapsed time is 5028.244784 seconds.
% Elapsed time is 5295.990748 seconds.
% Elapsed time is 5568.303231 seconds.
% Elapsed time is 5094.352764 seconds.
% Elapsed time is 5455.230147 seconds.
% Elapsed time is 5495.112906 seconds.
% Elapsed time is 5223.932153 seconds.
% Elapsed time is 5623.719860 seconds.
% Elapsed time is 5503.265266 seconds.
% Elapsed time is 5469.776985 seconds.
% Elapsed time is 5479.745458 seconds.
% Elapsed time is 5568.954163 seconds.
pool = gcp('nocreate');
delete(pool);

computeTimes = [4626.300634,4888.525983,5001.064181, 5430.624619,5522.823121,5629.031060,...
5686.836983,5697.228781,5775.514154,5777.833345,5818.694876 ,5856.220134,5487.011220,5251.232889,...
5497.733604,5694.623610,5649.478740,5566.005927,5571.943266,5705.395601,5776.357749,5729.238302,...
5690.881285, 5732.525444, 5028.244784, 5295.990748,5568.303231, 5094.352764, 5455.230147,...
5495.112906, 5223.932153, 5623.719860, 5503.265266, 5469.776985, 5479.745458,5568.954163];

mean(computeTimes)/60
median(computeTimes)/60

captureTimes     = ones(1,18);
nSamplesThinned = 20000;
nChains = 36;
distsL1a = zeros(nChains,nSamplesThinned);
permuts1 = zeros(nSamplesThinned,18);
refPos = 1:18;

for j = 1:nChains
orders1 = csvread(sprintf('Shalek13PreOrderedAV_Results_Orders_Chain%d.csv',j+300));
permStart1 = csvread(sprintf('Shalek13PreOrderedAV_Permutation%d.csv',j+300));
for k = 1:nSamplesThinned
    aa(permStart1(orders1(k,:))) = 1:18;
    permuts1(k,:) = aa;
    distsL1a(j,k) = pdist([aa; refPos],'cityblock');
end
end
csvwrite('Shalek13DistFromRefL1Clustg.csv',distsL1a);

set(0,'defaultAxesFontSize',6);
distsL1a = csvread('Shalek13DistFromRefL1Clustg.csv');
figure()
for j = 1:36
    ax(j) = subplot(6,6,j);
    histogram(distsL1a(j,10000:1:end),20,'BinLimits',[75 125])
    ylim([0,10000]);
end
set(gcf, 'PaperUnits', 'centimeters');
 x_width=17.6 ;y_width=17.6;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('Shalek13Hists','-dpdf','-r350')

%find corresponding PSMs
for j=1:36
    A = dlmread(sprintf('Shalek13PreOrderedAV_Results_Chain%d.csv',j+300),',',[9999 1 19999 nGenes]); 
    PSM_total = psm(A);
    csvwrite(sprintf('PSM_Shalek13g_%d.csv',j),PSM_total);
end

PSMs = zeros(142,142,36);
set(0,'defaultAxesFontSize',6);
for j=1:36
    A = csvread(sprintf('PSM_Shalek13g_%d.csv',j));
    PSMs(:,:,j) = A;
end
meanPSM = mean(PSMs,3);
[x,ind] = sort(sumClust(meanPSM,2:15));
save('Shalek13_PSMs.mat','PSMs');

figure()
for j=1:36
    ax(j) = subplot(6,6,j);
    imagesc(PSMs(ind,ind,j)); 
    colormap(gray);
    set(gca,'YDir','normal');
end
s2Pos = get(ax(36),'position');
s1Pos = get(ax(31),'position');
col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.01  s2Pos(2)  0.02  s1Pos(2)+6.9*s1Pos(3)]);
set(gcf, 'PaperUnits', 'centimeters');
 x_width=17.6 ;y_width=17.6;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('Shalek13PSMs','-dpdf','-r350')

alphas = [];
for j = 1:36
    File = sprintf('Shalek13PreOrderedAV_Results_Chain%d.csv',j+300);
    xx = dlmread(File,',',[0 0 20000-1 0]);
    alphas = [alphas xx];
end
csvwrite('alphasShalek13.csv',alphas);


%plot the summary PSM, you need to run Shalek13NonparametricSummaryPSM.R
set(0,'DefaultAxesFontSize',14);
load('npSummaryPSMShalek13.mat');
f=figure()
f.Units = 'centimeters';
imagesc(npSummaryPSMShalek13.weightedPSM(ind,ind));
set(gca,'YDir','normal')
colormap(gray);
col = colorbar;
col.FontSize= 14;
x_width=5 ;y_width=4;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
set(gcf,'PaperSize',[x_width y_width]);
print('PSM_PY_Shalek13','-dpdf','-r350')

%summary of orderings
captureTimes     = ones(1,18);
nSamplesThinned = 20000;
nChains = 36;
distsL1a = zeros(nChains,nSamplesThinned);
permuts1 = zeros(nSamplesThinned,18);
refPos = 1:18;

posits = [];
for j = 1:nChains
orders1 = csvread(sprintf('Shalek13PreOrderedAV_Results_Orders_Chain%d.csv',j+300));
permStart1 = csvread(sprintf('Shalek13PreOrderedAV_Permutation%d.csv',j+300));
for k = 10000:20:nSamplesThinned
    aa(permStart1(orders1(k,:))) = 1:18;
    posits = [posits;aa];
end
end
mdsSt= cmdscale(pdist(posits),2);
save('mdsSt_Shalek13.mat','mdsSt');
load('mdsSt_Shalek13.mat');
figure()
colo = ['r','c','m','g','y'];
[IDX, C, SUMD, D, MIDX] = kmedoids(mdsSt,4);
for k = 1:4
    k
scatter(mdsSt(IDX==k,1),mdsSt(IDX==k,2),100,'MarkerFaceColor',colo(k),...
'MarkerEdgeColor',colo(k),'MarkerEdgeAlpha',0.05,...
    'MarkerFaceAlpha',0.05);
hold on;
end

%compare orders obtained by GPseudoClust versus -Rank

positsRank = [];
for j = 1:3
orders1 = csvread(sprintf('Shalek13PreOrderedAV_Results_Orders_Chain%d.csv',j));
    permStart1 = csvread(sprintf('Shalek13PreOrderedAV_Permutation%d.csv',j));
for k = 10000:20:20000
    aa(permStart1(orders1(k,:))) = 1:18;
    positsRank = [positsRank;aa];
end
end
mdsStAll= cmdscale(pdist([posits(1:10:end,:);positsRank]),2);
save('mdsSt_Shalek13All.mat','mdsStAll');

load('mdsSt_Shalek13All.mat'); 
methodInd = [zeros(1,1804),ones(1,1503)];
figure() 
h0=scatter(mdsStAll(methodInd==0,1),mdsStAll(methodInd==0,2),100,'MarkerFaceColor','m',...
'MarkerEdgeColor','m','MarkerEdgeAlpha',0.05,...
    'MarkerFaceAlpha',0.05);
hold on;
h1=scatter(mdsStAll(methodInd==1,1),mdsStAll(methodInd==1,2),100,'MarkerFaceColor','b',...
'MarkerEdgeColor','b','MarkerEdgeAlpha',0.05,...
    'MarkerFaceAlpha',0.05);

legend([h0,h1],'GPseudoClust','GPseudoRank','Location','NorthOutSide');
x_width = 6;
y_width = 6;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
set(gcf,'PaperSize',[x_width y_width]);
print('Shalek13Comp','-dpdf','-r350')


computePositionsPseudoRank('Shalek13PreOrderedAV.csv',301:336,301:336,10000,...
    ones(1,18),'Shalek13Pos.mat');

computePositionsPseudoRank('Shalek13PreOrderedAV.csv',1:3,1:3,10000,...
    ones(1,18),'Shalek13PosRank.mat');

%add path to subaxis function
%https://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot
addpath(genpath('~/subaxis'));

clear all;
close all;
set(0,'DefaultAxesFontSize',16);
figure()
load('Shalek13Pos.mat');
ax(1) = subaxis(2,1,1,'SpacingHoriz', 0.025, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.1,...
         'MarginLeft',0.1,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(mean(posFreq,3));
set(gca,'Ydir','Normal');
xlabel('cell');
ylabel('position in order');
colormap(hot);
load('Shalek13PosRank.mat');
ax(2) = subaxis(2,1,2,'SpacingHoriz', 0.025, 'SpacingVert',0.05,'Padding', 0.01, 'MarginRight', 0.1,...
         'MarginLeft',0.1,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(mean(posFreq,3));
xlabel('cell');
ylabel('position in order');
set(gca,'Ydir','Normal');
colormap(hot);

s2Pos = get(ax(2),'position');
s1Pos = get(ax(2),'position');
x_width = 6;
y_width = 12;
col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.01  s2Pos(2)  0.05  s1Pos(2)+s1Pos(3)*0.86]);
col.FontSize=12;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
set(gcf,'PaperSize',[x_width y_width]);
print('Shalek13Comp2','-dpdf','-r350')



%turn around orders 
figure()
load('Shalek13Pos.mat');
ax(1) = subaxis(2,1,1,'SpacingHoriz', 0.025, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.1,...
         'MarginLeft',0.1,'MarginTop',0.1,'MarginBottom',0.1);
A = mean(posFreq,3);
Asymm = A;
Asymm(:,1:9) = (A(:,1:9)+A(:,18:-1:10));
Asymm(:,10:18) = (A(:,9:-1:1)+A(:,10:1:18));
ARev = zeros(size(A));
ARev(1:9,1:9) = Asymm(1:9,1:9);
ARev(10:18,10:18) = Asymm(10:18,10:18);
imagesc(ARev);
xlabel('cell');
ylabel('position in order');
set(gca,'Ydir','Normal');
colormap(hot);
load('Shalek13PosRank.mat');
ax(2) = subaxis(2,1,2,'SpacingHoriz', 0.025, 'SpacingVert',0.05,'Padding', 0.01, 'MarginRight', 0.1,...
         'MarginLeft',0.1,'MarginTop',0.1,'MarginBottom',0.1);

A = mean(posFreq,3);
Asymm = A;
Asymm(:,1:9) = (A(:,1:9)+A(:,18:-1:10));
Asymm(:,10:18) = (A(:,9:-1:1)+A(:,10:1:18));
ARev = zeros(size(A));
ARev(1:9,1:9) = Asymm(1:9,1:9);
ARev(10:18,10:18) = Asymm(10:18,10:18);
imagesc(ARev);
set(gca,'Ydir','Normal');
colormap(hot);
xlabel('cell');
ylabel('position in order');

s2Pos = get(ax(2),'position');
s1Pos = get(ax(2),'position');
x_width = 6;
y_width = 12;
col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.01  s2Pos(2)  0.05  s1Pos(2)+s1Pos(3)*0.86]);
col.FontSize=12;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
set(gcf,'PaperSize',[x_width y_width]);
print('Shalek13ComNonSymm','-dpdf','-r350')


%reversed, log-scale
figure()
load('Shalek13Pos.mat');
ax(1) = subaxis(2,1,1,'SpacingHoriz', 0.025, 'SpacingVert',0.1,'Padding', 0.01, 'MarginRight', 0.1,...
         'MarginLeft',0.1,'MarginTop',0.1,'MarginBottom',0.1);
A = mean(posFreq,3);
Asymm = A;
Asymm(:,1:9) = (A(:,1:9)+A(:,18:-1:10))*0.5;
Asymm(:,10:18) = (A(:,9:-1:1)+A(:,10:1:18))*0.5;
ARev = zeros(size(A));
ARev(1:9,1:9) = Asymm(1:9,1:9);
ARev(10:18,10:18) = Asymm(10:18,10:18);
imagesc(log2(ARev));
xlabel('cell');
ylabel('position in order');
set(gca,'Ydir','Normal');
colormap(hot);
load('Shalek13PosRank.mat');
ax(2) = subaxis(2,1,2,'SpacingHoriz', 0.025, 'SpacingVert',0.05,'Padding', 0.01, 'MarginRight', 0.1,...
         'MarginLeft',0.1,'MarginTop',0.1,'MarginBottom',0.1);

A = mean(posFreq,3);
Asymm = A;
Asymm(:,1:9) = (A(:,1:9)+A(:,18:-1:10))*0.5;
Asymm(:,10:18) = (A(:,9:-1:1)+A(:,10:1:18))*0.5;
ARev = zeros(size(A));
ARev(1:9,1:9) = Asymm(1:9,1:9);
ARev(10:18,10:18) = Asymm(10:18,10:18);
imagesc(log2(ARev));
xlabel('cell');
ylabel('position in order');
set(gca,'Ydir','Normal');
colormap(hot);

s2Pos = get(ax(2),'position');
s1Pos = get(ax(2),'position');
x_width = 6;
y_width = 12;
col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.01  s2Pos(2)  0.05  s1Pos(2)+s1Pos(3)*0.86]);
col.FontSize=12;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
set(gcf,'PaperSize',[x_width y_width]);

print('Shalek13ComNonSymmLog','-dpdf','-r350')
