function[postSim] = psm(clusterAllocations)
% compute the posterior similarity matrix for cluster allocation samples 
% from the posterior distribution
%input: n_samples x n matrix, each row is one sample of cluster allocations
%after burn-out
[M,N] = size(clusterAllocations);
postSim = zeros(N,N);
for i = 1:N
    for j = 1:(i-1)
        for m = 1:M
            if clusterAllocations(m,i) == clusterAllocations(m,j)
                postSim(i,j) = postSim(i,j)+1/M;
            end
        end
    end
    postSim(i,i)=0.5;
end
postSim = postSim+postSim';
end