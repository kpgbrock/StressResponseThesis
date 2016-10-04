function [results] = pairSwapMutVHL()

load vhl;
addpath('..');
x = 63:204;
seq = vhl.orgSeq(x);

[copt,V] = pholder(vhl.mutSeq{19}(x));
temp19 = V*copt;
[tempCrys,~] = createCrystalBurialTraceChain('1VCB',seq,'C');

%aa = {'I','V','L','F','C','M','A','G','T','W','S','Y','P','H','E','Q','D','N','K','R'};

nRes = length(x);

N = nRes*(nRes-1)/2;
results.pairSwap = zeros(N,2);
results.bt = zeros(nRes,N);
results.corrCrystal = zeros(N,1);
results.corr19 = zeros(N,1);

z = 1;
for i=1:nRes
    disp(i);
    for j=(i+1):nRes
        seqNew = seq;
        seqNew([i,j]) = seqNew([j,i]);
        [copt,V] = pholder(seqNew);
        
        results.pairSwap(z,:) = [i,j];
        results.bt(:,z) = V*copt;
        results.corrCrystal(z) = corr(results.bt(:,z),tempCrys);
        results.corr19(z) = corr(results.bt(:,z),temp19);
        z = z+1;
    end
end
end