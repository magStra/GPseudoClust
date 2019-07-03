addpath(genpath('../'));
addpath('~/subaxis');
%example of a data set without capture times
CT     = [repmat(1,1,35)];
fHandle          = @GPseudoClust2;  
fileName         = 'mesCC.csv';
nSamples         = 4000; % Number of thinned samples
verbose          = false; % Whether or not to print output to screen
inputSeed        = NaN;     % Set to NaN if you do not want a seed, this will allocate a seed
permuteData      = true;%if true, the sampler starts from a randomly permuted order
b                = 0.01;
adjustForCellSize = false;
parpool(4);
parfor j = 1:4
    tic
    feval(fHandle, fileName, j, nSamples, verbose,inputSeed,permuteData,CT,...
                b,adjustForCellSize,1:35);
    toc
%     Elapsed time is 5374.177827 seconds.
%     Elapsed time is 5569.817745 seconds.
%     Elapsed time is 5593.426180 seconds.
%     Elapsed time is 5687.290024 seconds.
     A = dlmread(sprintf('mesCC_Results_Chain%d.csv',j),',',[1999 1 3999 600]); 
    PSM_total = psm(A);
        csvwrite(sprintf('PSM_mes%d.csv',j),PSM_total);
end
computeTimes = [5374.177827,5569.817745,5593.426180,5687.290024];
median(computeTimes)/60
A = zeros(600,600,3);
for j = 1:3
    A(:,:,j)=csvread(sprintf('PSM_mes%d.csv',j));
end
%ordering of the rows and columns of the matrices to make the cluster
%structure visible. Please note that the labels in the summary clustering
%may swap across runs, and therefore the blocks of genes with high
%co-clustering probabilities may be permuted, and the plot therefore look
%slightly different for each run (blocks differently arranged), even if the PSMs themselves are almost
%identical. 
set(0,'defaultAxesFontSize',8);
Amean = mean(A,3);
cl = sumClust(Amean,2:12);
[~,ind] = sort(cl);

for j=1:3
    ax(j) = subaxis(1,3,j,'SpacingHoriz', 0.03, 'SpacingVert',0.01,'Padding', 0.01, 'MarginRight', 0.06,...
         'MarginLeft',0.05,'MarginTop',0.05,'MarginBottom',0.1);
    imagesc(A(ind,ind,j));
    set(gca,'YDir','normal')
    colormap(gray);
   
end
s2Pos = get(ax(3),'position');
s1Pos = get(ax(3),'position');
col = colorbar('Position', [s2Pos(1)+s2Pos(3)+0.01  s2Pos(2)  0.02  s1Pos(2)+s1Pos(3)*2.83])

set(gcf, 'PaperUnits', 'centimeters');
 x_width=17.6 ;y_width=4.5;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('PSMS_mesNoSubs','-dpdf','-r350')


