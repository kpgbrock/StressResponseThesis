function [results] = pHRandomControl()

addpath('..');
load yCytoExp;
yc = yCytoExp;
N = length(yc.id);
pHRange = 1:0.1:14;

[~,residuesPos,residuesNeg] = findpKaWiki('');
res = strcat(residuesPos,residuesNeg);
resCounts = zeros(length(res),1);

results.seq = cell(N,1);
results.pHRange = pHRange;
results.charge = zeros(length(pHRange),N);
results.pI = zeros(N,1);
results.info = 'Randomized distribution of charged residues that matches yeast proteome.';

for i=1:N
    for j=1:length(resCounts)
        s = yc.sequence{i};
        resCounts(j) = resCounts(j) + length(s(s==res(j)));
        
    end
end

results.resCounts = resCounts;

newSeq = cell(N,1);
for i=1:length(resCounts)
    for j=1:resCounts(i)
       z = randi(N);
       newSeq{z} = strcat(newSeq{z},res(i));
    end
end

results.seq = newSeq;

for i=1:N
    seq = newSeq{i};
    [fcharge,pI] = findPICurveWikiNoEnds(seq,pHRange);
    results.charge(:,i) = fcharge;
    results.pI(i) = pI;
end

end