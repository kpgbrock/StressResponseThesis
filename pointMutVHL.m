function [results] = pointMutVHL()

load vhl;
addpath('..');

x = 63:204;
seq = vhl.orgSeq(x);

[copt,V] = pholder(vhl.mutSeq{19}(x));
temp19 = V*copt;
[tempCrys,~] = createCrystalBurialTraceChain('1VCB',seq,'C');

aa = {'I','V','L','F','C','M','A','G','T','W','S','Y','P','H','E','Q','D','N','K','R'};

N = length(x);
results.aa = aa;
results.ptMut = zeros(20*N,2);
results.bt = zeros(N,20*N);
results.corrCrystal = zeros(20*N,1);
results.corr19 = zeros(20*N,1);

for i=1:N
    disp(i);
    seqNew = seq;
    for j=1:20
        
        seqNew(i) = aa{j};
        [copt,V] = pholder(seqNew);
        z = (i-1)*20 + j;
        results.ptMut(z,:) = [i,j];
        results.bt(:,z) = V*copt;
        results.corrCrystal(z) = corr(results.bt(:,z),tempCrys);
        results.corr19(z) = corr(results.bt(:,z),temp19);
        
    end
end
end