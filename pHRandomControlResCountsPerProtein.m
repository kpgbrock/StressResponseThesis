function [results] = pHRandomControlResCountsPerProtein()

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

resL = zeros(N,1);
resChar = '';

for i=1:N
    s = yc.sequence{i};
    for j=1:length(resCounts)  
        nR = length(s(s==res(j)));
        resL(i) = resL(i) + nR;
        resChar = strcat(resChar,s(s==res(j)));
    end
end

resChar = resChar(randperm(length(resChar)));
newSeq = cell(N,1);
for i=1:N
    newSeq{i} = resChar(1:resL(i));
    
    resChar = resChar((resL(i)+1):end);
    %disp(newSeq{i});
    [fcharge,pI] = findPICurveWikiNoEnds(newSeq{i},pHRange);
    results.charge(:,i) = fcharge;
    if isempty(pI)
        pI = NaN;
    end
    results.pI(i) = pI;
end
results.seq = newSeq;

end