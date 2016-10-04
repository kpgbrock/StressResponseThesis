function [results] = isingModelIntrxns(pairIDs, phNeutral,phAcidic, phycfull, fxEnergy)

%load phycfull;

% If no energy function specified, use default
% If we want to use yeast pH/charge data, load phycfull; else, we passed in
% the data structure as the fifth argument
if nargin == 3
    load phycfull;
    fxEnergy = @getEnergyQL;
elseif nargin == 4
    fxEnergy = @getEnergyQL;
end


[numPairs,~] = size(pairIDs);
pChargeNeutral = zeros(numPairs,2);
pChargeAcidic = zeros(numPairs,2);
pLength = zeros(numPairs,2);

neutID = find(phycfull.pHRange == phNeutral);
acidID = find(phycfull.pHRange == phAcidic);

results.eTotalNeutral = 0;
results.eTotalAcidic = 0;
results.dE = 0;
results.eComponentsNeutral = zeros(numPairs,1);
results.eComponentsAcidic = zeros(numPairs,1);

for i=1:numPairs
    
    pChargeNeutral(i,1) = phycfull.charge(neutID,pairIDs(i,1));
    
    pChargeNeutral(i,2) = phycfull.charge(neutID,pairIDs(i,2));
    pChargeAcidic(i,1) = phycfull.charge(acidID,pairIDs(i,1));
    pChargeAcidic(i,2) = phycfull.charge(acidID,pairIDs(i,2));
    
    pLength(i,1) = phycfull.seqLength(pairIDs(i,1));
    pLength(i,2) = phycfull.seqLength(pairIDs(i,2));
    
    
    
    

end

results.eComponentsNeutral = fxEnergy(pChargeNeutral, pLength);
results.eComponentsAcidic = fxEnergy(pChargeAcidic, pLength);
results.eTotalNeutral = sum(results.eComponentsNeutral);
results.eTotalAcidic = sum(results.eComponentsAcidic);
results.dE = results.eTotalAcidic - results.eTotalNeutral;

end

function [energyvals] = getEnergyQL(qvals, pLength)

z = qvals./pLength;
[r,~] = size(z);
energyvals = zeros(r,1);
for i=1:r
    energyvals(i) = z(i,1)*z(i,2);
end

end