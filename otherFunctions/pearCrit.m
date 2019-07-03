function outcome = pearCrit(cl,PSM)
%PEAR criterion by Fritsch and Ickstadt (2009). Bayesian analysis. 
%Improved criteria for clustering based on the posterior similarity matrix. 
if size(cl,1) > 1
    cl = cl';
end
mat = psm(cl).*PSM;
sumIpi = sum(sum(tril(mat,-1)));
n = size(PSM,1);
no2 = nchoosek(n,2);
sumpij = sum(sum(tril(PSM,-1)));
sumIij = sum(sum(tril(psm(cl),-1)));
correc = (sumIij*sumpij)/no2;
outcome = (sumIpi - correc)/ (0.5*(sumpij+sumIij)-correc);
end
