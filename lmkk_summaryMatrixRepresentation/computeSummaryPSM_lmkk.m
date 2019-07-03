function [clusterSolution ncl] = computeSummaryPSM_lmkk(PSMs,nc)
%This computes summary matrix representations from multiple PSMs using
%localised multiple kernel learning. 

%Gönen, M. and Margolin, A.A. (2014). Localized data fusion for kernel k-means
%clustering with application to cancer biology. In Advances in Neural Information
%Processing Systems 27, pages 1305-1313.

%This function requires the lmkk code accompanying the above publication and 
%available at https://github.com/mehmetgonen/lmkkmeans

%and the Mosek optimisation software (https://www.mosek.com/)

%The number of clusters is estimated before using the the criterion used in
%SIMLR (see estimateNumberClustersSIMLR.m)

%input: PSMs: 3-dimensional array of size nSamples x nSamples x nViews
% ncl: range of cluster numbers to be considered

addpath('~/lkkmeans')%add path to folder containing lmkk code
addpath(genpath('~/mosek'));%add path to folder containing the mosek software
meanPSM = mean(PSMs,3);

if length(nc)>1
    [K1 K2] = estimateNumberClustersSIMLR_GPseudoClust_lmkk(meanPSM,2:10);
    [KK2 ind2] = max(K2);
    ncl = nc(ind2); 
else
    ncl=nc;
end
parameters = struct();
parameters.iteration_count = 10;
d = size(PSMs,1);
PSM = zeros(d,d);
parameters.cluster_count = ncl;
clusterSolution = lmkkmeans_train(PSMs, parameters);
weights = clusterSolution.Theta;
K = clusterSolution.K_Theta;
for j = 1:d
    for i = 1:d
        PSM(i,j) = K(i,j)/sum(weights(i,:).*weights(j,:));
    end
end
clusterSolution.PSM = PSM;
disp(ncl)
end
