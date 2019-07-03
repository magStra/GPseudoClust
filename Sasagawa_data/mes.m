addpath(genpath('../'));
CT     = [repmat(1,1,35)];
fHandle          = @GPseudoClust2;  
fileName         = 'mesCC.csv';
nSamples         = 4000; % Number of thinned samples
verbose          = false; % Whether or not to print output to screen
inputSeed        = NaN;     % Set to NaN if you do not want a seed, this will allocate a seed
permuteData      = true;%if true, the sampler starts from a randomly permuted order
b                = 0.01;
adjustForCellSize = false;
PSMS = zeros(600,600,36);
tic
parpool(12);
%adjust to number of cores
parfor j = 1:36
 subS = randsample(1:35,15);
 tic
  feval(fHandle, fileName, j+100, nSamples, verbose,inputSeed,permuteData,ones(1,15),...
               b,adjustForCellSize,subS);
    A = dlmread(sprintf('mesCC_Results_Chain%d.csv',j+100),',',[1999 1 3999 600]);
    toc
%  Elapsed time is 3872.693006 seconds.
% Elapsed time is 3915.107841 seconds.
% Elapsed time is 4010.744007 seconds.
% Elapsed time is 4046.912471 seconds.
% Elapsed time is 4081.392788 seconds.
% Elapsed time is 4099.033860 seconds.
% Elapsed time is 4166.531643 seconds.
% Elapsed time is 4343.430221 seconds.
% Elapsed time is 4363.376893 seconds.
% Elapsed time is 4373.851505 seconds.
% Elapsed time is 4735.719973 seconds.
% Elapsed time is 4770.471853 seconds.
% Elapsed time is 3586.192886 seconds.
% Elapsed time is 3720.569925 seconds.
% Elapsed time is 3803.804025 seconds.
% Elapsed time is 4062.283697 seconds.
% Elapsed time is 4094.226311 seconds.
% Elapsed time is 4076.786074 seconds.
% Elapsed time is 4346.329041 seconds.
% Elapsed time is 3707.383364 seconds.
% Elapsed time is 4170.610336 seconds.
% Elapsed time is 4208.394785 seconds.
% Elapsed time is 4710.512781 seconds.
% Elapsed time is 4383.019709 seconds.
% Elapsed time is 4024.159621 seconds.
% Elapsed time is 4389.684194 seconds.
% Elapsed time is 4335.734571 seconds.
% Elapsed time is 4230.236594 seconds.
% Elapsed time is 4347.385779 seconds.
% Elapsed time is 4011.637565 seconds.
% Elapsed time is 4070.901405 seconds.
% Elapsed time is 3987.546947 seconds.
% Elapsed time is 4144.279008 seconds.
% Elapsed time is 3672.056567 seconds.
% Elapsed time is 4255.145487 seconds.
% Elapsed time is 4374.611701 seconds.
% Elapsed time is 13293.807802 seconds.


    PSMS(:,:,j) = psm(A);
end
toc
save('PSMs_mes_sub.rda','PSMS');
computeTimes = [3872.693006,3915.107841, 4010.744007,...
4046.912471, 4081.392788,4099.033860,4166.531643,4343.430221,4363.376893,...
4373.851505,4735.719973,4770.471853,3586.192886,3720.569925,3803.804025,4062.283697,4094.226311,... 
4076.786074,4346.329041,3707.383364,4170.610336,4208.394785,4710.512781,4383.019709,... 
4024.159621,4389.684194,4335.734571,4230.236594,4347.385779,4011.637565,4070.901405,... 
3987.546947,4144.279008,3672.056567,4255.145487,4374.611701,13293.807802];

median(computeTimes)/60

alphasSasagawa = zeros(2000,36);
for j = 1:36
    alphasSasagawa(:,j) = dlmread(sprintf('mesCC_Results_Chain%d.csv',j+100),',',[0 0 1999 0]);
end
csvwrite('alphasSasagawa.csv',alphasSasagawa);

%plotting
load('PSMs_mes_sub.mat');
clusterSolution = computeSummaryPSM_lmkk(PSMS,2:12);
figure()
ax(4) = subaxis(1,4,4, 'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(clusterSolution.PSM(ind,ind));
set(gca,'YDir','normal')
title('lmkk');
colormap(gray);
ax(3) = subaxis(1,4,3, 'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
A = mean(PSMS,3);
imagesc(A(ind,ind));
set(gca,'YDir','normal')
colormap(gray);
title('mean PSM');
%processing the results from PY
load('npSummaryPSMSasagawa.mat');
ax(1) = subaxis(1,4,1,'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(npSummaryPSMSasagawa.weightedPSM(ind,ind));
set(gca,'YDir','normal')
title('PY and PEAR');
ax(2) = subaxis(1,4,2,'SpacingHoriz', 0.025, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.05,...
         'MarginLeft',0.02,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(npSummaryPSMSasagawa.weightedPSM_DP(ind,ind));
set(gca,'YDir','normal')
title('DPM and PEAR');
s2Pos = get(ax(4),'position');
s1Pos = get(ax(4),'position');
col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.01  s2Pos(2)  0.02  s1Pos(2)+s1Pos(3)*3.5])

set(gcf, 'PaperUnits', 'centimeters');
 x_width=17.6 ;y_width=4;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('PSMS_mes','-dpdf','-r350')





