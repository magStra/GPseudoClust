addpath(genpath('../'));
set(0,'defaultFigureUnits','centimeters');
%sim2Large(sprintf('simLarge'));
fHandle          =  @GPseudoClust2;  
fileName         = 'simLarge.csv';
nSamples         = 5000;  
verbose          = false; 
inputSeed        = NaN;    
permuteData      = true;
b = 0.01;
adjustForCellSize = false;
A = csvread(fileName);
[nGenes nCells] = size(A);

captureTimes     = [ones([1,30]),repmat(2,[1,30]),repmat(3,[1,30])];
a = nCells/3;
PSMs = zeros(nGenes,nGenes,96);
parpool(12);
parfor jj = 1:96
 ab = [randsample(1:a,30),randsample((a+1):(2*a),30),randsample((2*a+1):(3*a),30)];
 feval(fHandle, fileName, jj+300, nSamples, verbose,inputSeed,...
                permuteData,captureTimes,b,adjustForCellSize,ab);
PSMs(:,:,jj) = psm(dlmread(sprintf('simLarge_Results_Chain%d.csv',jj+300),',',[2449 1 4999 nGenes]));
end
pool = gcp('nocreate');
delete(pool);
save('PSMsLarge30.mat','PSMs');

