function [ppath tau X] = ProposingPseudoTimes(inputDataStructure,regular,n0,n3,n3a)
%This function samples from the proposal distribution of the orders
%and returns the corresponding order and  pseudotimes.
%This function is taken from the GPseudoRank code. 
%https://github.com/magStra/GPseudoRank
%Strauss, M. et al. (2018). GPseudoRank: a permutation sampler for single cell
%orderings. Bioinformatics, page bty664.

X           = sum(inputDataStructure.pp<=rand);
nFeatures   = inputDataStructure.nFeatures;
order       = inputDataStructure.order;
switch X
    case 0
        % swaps of neighbours
        ppath       = inputDataStructure.order;
        aaa         = max(n0,1);
        bbb         = randi(aaa);
        ii1         = randsample(nFeatures-1,bbb);
        for j = 1:bbb
            i1 = ii1(j);
            i2 = i1+1;
            helpVar     = ppath([i1,i2]);
            ppath(i1)   = helpVar(2);
            ppath(i2)   = helpVar(1);
        end
   
    case 1
         %second case, swaps of similar elements
         Y         = sum(inputDataStructure.cumMoveProb1 <= rand)+1;
         ppath     = inputDataStructure.order;
         m        = mod(Y,nFeatures);
         n         = ceil(Y/nFeatures);
         if m == 0
             m = nFeatures;
         end
         i1         = find(ppath == n);
         i2         = find(ppath ==  m);
         ppath(i1)  = order(i2);
         ppath(i2)  = order(i1);
       
      case 2% third case, reversals
         Y         = sum(inputDataStructure.cumMoveProb2 <= rand)+1;
         ppath      = inputDataStructure.order;
         m          = mod(Y,nFeatures);
         n          = ceil(Y/nFeatures);
         if m == 0
             m = nFeatures;
             n = n-1;
         end
         i1         = find(ppath == n);
         i2         = find(ppath ==  m);
         ppath(min(i1,i2):max(i1,i2))  = inputDataStructure.order(max(i1,i2):-1:min(i1,i2)); 
     
         
        case 3%random permutations
        ppath     = inputDataStructure.order;
        nPerm     = randi(max(1,n3));
        %for each such permutation draw the length
        % we use a length of at least 3
        lPerms  = randi([3 max(3,n3a)],[1 nPerm]);
        for i = 1:nPerm
        %now for each permutation pick the first index
            i1 = randi(nFeatures-lPerms(i)+1);
            i2 = i1 + lPerms(i)-1;
            xx = ppath(i1:i2);
            helpVar = xx(randperm(i2-i1+1));
            ppath(i2:-1:i1) = helpVar;
        end 
end 
%now convert path to pseudo-time
dataProposedOrder = inputDataStructure.data(:,ppath);
if regular == false
    tau1 = sqrt(sum((diff(dataProposedOrder,1,2).^2),1));
    tau = cumsum([0 tau1])/sum(tau1);
else
    tau = (0.5:1:(nFeatures-0.5))/nFeatures;
end
if any(unique(inputDataStructure.order) ~= unique(ppath))
    'error -- new order is not a permutation of old one'
end 
if any(isnan(ppath))
   'error: NaN in order'
end
end

