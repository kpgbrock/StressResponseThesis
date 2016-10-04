function [results] = phRandSampleBothNeutralY2H(nTrials)

load y2hdata;
load phycfull;

setOfIntrxns = nchoosek(y2hdata.idxIndProtNS,2);
numPossInt = length(setOfIntrxns);
N = length(y2hdata.idxIndProtNS);

results.numOn = zeros(nTrials,1);
results.y2hChLen = 0;

z5 = find(phycfull.pHRange == 5);
z7 = find(phycfull.pHRange == 7);

for i=1:nTrials
    z = randperm(numPossInt);
    %z = z(1:length(y2hdata.idxCytPairs));
    z = z(1:length(y2hdata.idxNonSelf));
%     if mod(i,500) == 0
%         disp(i);
%     end
    lInt = setOfIntrxns(z,:);
    results.numOn(i) = length(find((abs(phycfull.charge(z5,lInt(:,1))) < abs(phycfull.charge(z7,lInt(:,1)))) & (abs(phycfull.charge(z5,lInt(:,2))) < abs(phycfull.charge(z7,lInt(:,2))))));
end

end