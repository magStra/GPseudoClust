%
clear all
close all
addpath(genpath('../'));
%add path to SIMLR
addpath(genpath('~/SIMLR-SIMLR'));

%Simulated data set 1
data = csvread('simDataClust.csv');
true_labs =  [repmat(1,1,8),repmat(2,1,4),repmat(3,1,12),repmat(4,1,16),repmat(5,1,12)];
%Estimate number of clusters
[K1, K2] = Estimate_Number_of_Clusters_SIMLR(data,2:10);
[x,k1] = min(K1);
[x,k2] = min(K2);
k1
k2
%both are 4, so it is 5 clusters
[y, S, F, ydata,alpha] = SIMLR(data,5);
NMI_simDAta1 = Cal_NMI(y,true_labs);
%save cluster allocations
csvwrite('sim1_SIMLR.csv',y);

%Simulated data set 2
data = csvread('simDataClust2.csv');
true_labs =  [repmat(1,1,8),repmat(2,1,4),repmat(3,1,12),repmat(4,1,16),repmat(5,1,12)];
%Estimate number of clusters
[K1, K2] = Estimate_Number_of_Clusters_SIMLR(data,2:10);
[x,k1] = min(K1);
[x,k2] = min(K2);
k1
k2
%both are 6, so it is 7 clusters
[y, S, F, ydata,alpha] = SIMLR(data,7);
NMI_simData2 = Cal_NMI(y,true_labs);
%save down cluster allocations
csvwrite('sim2_SIMLR.csv',y);


%Shalek data
addpath('../../preprocessedScRNASeq/');
Data = importdata('Shalek.csv');
data = Data.data(:,125:307);
[K1, K2] = Estimate_Number_of_Clusters_SIMLR(data,2:15);
[x,k1] = min(K1);
[x,k2] = min(K2);
k1
k2
%both are 11, so it is 12 clusters
[y, S, F, ydata,alpha] = SIMLR(data,12);
%save down cluster allocations
csvwrite('Shalek_SIMLR.csv',y);



