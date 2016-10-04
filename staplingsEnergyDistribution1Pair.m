function [results] = staplingsEnergyDistribution1Pair(seqA,seqB,numStaples,numTrials)

nA = length(seqA);
nB = length(seqB);

results.seqA = seqA;
results.seqB = seqB;
results.numStaples = numStaples;
results.numTrials = numTrials;
results.emin = zeros(numTrials,1);


for i=1:numTrials
    disp(i);
    pairs = zeros(numStaples,2);
    for j=1:numStaples
       pairs(j,1) = randi(nA);
       pairs(j,2) = randi(nB);
    end
    [~,~,~,emin,~,~,~,~] = pholderStapleAB(seqA,seqB,pairs);
    results.emin(i) = emin;
end
end