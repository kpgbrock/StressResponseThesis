function [results] = growMultiPolymer()

nC = 50;
nProt = 2;
cutOffLength = 10;
N = 100;

numConf = nC*ones(nProt,1);

results.AA = zeros(N*length(numConf)*2,1);
results.AB = zeros(N*length(numConf)*2,1);

for i=1:N
    [polymer] = calcPolymerizationTime(numConf,cutOffLength);
    results.AA(((i-1)*(nC*nProt)+1):(i*(nC*nProt))) = polymer.selfself(:,end);
    results.AB(((i-1)*(nC*nProt)+1):(i*(nC*nProt))) = polymer.AB(:,end);
end
end