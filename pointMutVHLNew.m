function [results] = pointMutVHLNew()

load vhlUP;
addpath('..');

x = 63:204;
seq = vhlUP.orgSeq(x);

[copt,V] = pholder(vhlUP.mutSeq{19}(x));
temp19 = V*copt;
[tempCrys,~] = createCrystalBurialTraceChain('1VCB',seq,'C');

aa = int2aa(1:20);

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
        
        seqNew(i) = int2aa(j);
        [copt,V] = pholder(seqNew);
        z = (i-1)*20 + j;
        results.ptMut(z,:) = [i,j];
        results.bt(:,z) = V*copt;
        results.corrCrystal(z) = corr(results.bt(:,z),tempCrys);
        results.corr19(z) = corr(results.bt(:,z),temp19);
        
    end
end
end