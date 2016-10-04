function [eigVals,bModes] = eigenmodeAnalysisVHL(x)

load vhl;
addpath('..');

seqs = cell(21,1);
seqs{1} = vhl.orgSeq(x);
for i=1:20
    seqs{i+1} = vhl.mutSeq{i}(x);
end

eigVals = zeros(length(x),21);
bModes = cell(21,1);
for i=1:21
    [eps,~,psi] = pholderEig(seqs{i});
    eigVals(:,i) = eps;
    bModes{i} = psi;
end
end