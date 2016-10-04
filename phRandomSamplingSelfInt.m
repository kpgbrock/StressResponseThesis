function [results] = phRandomSamplingSelfInt(nTrials)

load y2hdata;
load phycfull;


nSI = length(y2hdata.idxSelfInt);
N = length(y2hdata.idxCytPairs);

z5 = find(phycfull.pHRange == 5);
z7 = find(phycfull.pHRange == 7);

results.sum57 = zeros(nSI,nTrials);
results.absSum57 = zeros(nSI,nTrials);
results.pINorm = zeros(nSI,nTrials);
results.pI = zeros(nSI,nTrials);

for i=1:nTrials
   z = randperm(N);
   z=z(1:nSI);
   
   results.sum57(:,i) = (phycfull.charge(z5,z)+phycfull.charge(z7,z))'./phycfull.seqLength(z);
   results.absSum57(:,i) = (abs(phycfull.charge(z5,z))+abs(phycfull.charge(z7,z)))'./phycfull.seqLength(z);
   results.pI(:,i) = (phycfull.pI(z));
   results.pINorm(:,i) = phycfull.pI(z)./phycfull.seqLength(z);
end

end