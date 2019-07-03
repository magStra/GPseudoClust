function sumClust = sumClust(PSM,nc)
%psm...posterior similarity matrix
%nc...cluster numbers to be tested
pears = zeros(length(nc),1);
sumClusts = zeros(length(nc),size(PSM,1));
Z=linkage(1-PSM,'complete');
for j = 1:length(nc)
    sumClusts(j,:)= cluster(Z,'maxclust',nc(j));
    pears(j) = pearCrit(sumClusts(j,:),PSM);
end
[x ind] = max(pears);
sumClust = sumClusts(j,:);
end

