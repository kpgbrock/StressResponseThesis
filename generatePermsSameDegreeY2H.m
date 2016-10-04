function [results] = generatePermsSameDegreeY2H(numPerm)

load y2hdata;
load phycfull;

binZ = y2hdata.idxNonSelf;

actResults = isingModelIntrxns(binZ,7,5,phycfull);

results.actdE = actResults.dE;
results.actETotalNeutral = actResults.eTotalNeutral;
results.actETotalAcidic = actResults.eTotalAcidic;
results.actEComponentsNeutral = actResults.eComponentsNeutral;
results.actEComponentsAcidic = actResults.eComponentsAcidic;

for i=1:length(binZ)
    if binZ(i,1) > binZ(i,2)
        binZ(i,:) = [binZ(i,2) binZ(i,1)];
    end
end
binZ = sortrows(binZ,1);
indZ = y2hdata.idxIndProtNS;
N = length(indZ);
[NInt,~] = size(binZ);

for i=1:numPerm
    newProts = y2hdata.idxIndProtNS(randperm(N));
    binZNew = binZ;
    
    for j=1:NInt
        binZNew(j,1) = newProts(indZ==binZNew(j,1));
        binZNew(j,2) = newProts(indZ==binZNew(j,2));
    end
    
    tempResults = isingModelIntrxns(binZNew,7,5,phycfull);
    
    results.dEHist(i) = tempResults.dE;
    results.eTotalNeutral(i) = tempResults.eTotalNeutral;
    results.eTotalAcidic(i) = tempResults.eTotalAcidic;
    results.eCompNeutral{i} = tempResults.eComponentsNeutral;
    results.eCompAcidic{i} = tempResults.eComponentsAcidic;
end
        


end

