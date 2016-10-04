function [results] = pairSwapMutVHLUntruncated()

load vhlUP;
%addpath('..');
x = 63:204;
seq = vhlUP.orgSeq;

nRes = length(vhlUP.orgSeq);

N = nRes*(nRes-1)/2;
results.pairSwap = zeros(N,2);
results.bt = zeros(length(x),N);

z = 1;
for i=1:nRes
    disp(i);
    for j=(i+1):nRes
        seqNew = seq;
        seqNew([i,j]) = seqNew([j,i]);
        [copt,V] = pholder(seqNew(x));
        
        results.pairSwap(z,:) = [i,j];
        results.bt(:,z) = V*copt;
       
        z = z+1;
    end
end
end