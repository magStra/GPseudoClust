function [co,noCo] = coClustGeneLists(weightedPSM,genes,pMin,pMax)
%finds groups of genes co-clustered with a probability of at least pMin,
% and groups of genes co-clustered with a probability of less than pMax
% other inputs: weighted PSM (computed for instance using the 'PY+PEAR' method)
% genes: cell array of gene names
weightedPSMa =  weightedPSM > pMin;
weightedPSMb =  weightedPSM < pMax;
ng = size(weightedPSM,1);
neighbours = {};
notNeighbours = {};
k = 0;
k1 = 0;
for j = 1:ng
    a = unique(neighbors(graph(weightedPSMa),j));
    b = unique(neighbors(graph(weightedPSMb),j));
    if length(a)>= 2
        k = k+1;
        neighbours{k} =a;
    end
    if length(b)>= 2
        k1 = k1+1;
        notNeighbours{k1} =b;
    end
end
%remove duplicates
helpCell = cellfun(@(x) num2str(x(:)'),neighbours,'UniformOutput',false);
[~,ind,~] = unique(helpCell);
neighbours = neighbours(ind);
helpCell = cellfun(@(x) num2str(x(:)'),notNeighbours,'UniformOutput',false);
[~,ind,~] = unique(helpCell);
notNeighbours = notNeighbours(ind);
%checking if any of those is a subset of any other
indicator = true(1,length(neighbours));
for j =1:length(neighbours)
    a = neighbours{j};
    for k = 1:length(neighbours)
        b = neighbours{k};
        if j ~= k && all(ismember(b,a))
            indicator(k) = false;
        end
    end            
end

neighbours = neighbours(indicator);

indicator = true(1,length(notNeighbours));
for j =1:length(notNeighbours)
    a = notNeighbours{j};
    for k = 1:length(notNeighbours)
        b = notNeighbours{k};
        if j ~= k && all(ismember(b,a))
            indicator(k) = false;
        end
    end            
end
notNeighbours = notNeighbours(indicator);

co = {};
ij = 0;
for j = 1:length(neighbours)
    if length(neighbours{j}) > 1
        A = {};
        ij = ij +1;
        for k = 1:length(neighbours{j})
            A{k} = genes{neighbours{j}(k)};
        end
        co{ij} = A;
    end
end
noCo = {};
ij = 0;
for j = 1:length(notNeighbours)
    if length(notNeighbours{j}) > 1
        A = {};
        ij = ij +1;
        for k = 1:length(notNeighbours{j})
            A{k} = genes{notNeighbours{j}(k)};
        end
        noCo{ij} = A;
    end
end

end

