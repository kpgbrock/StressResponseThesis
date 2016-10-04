function [results] = stapleExp4112013()

load yCyto100300;
seqA = yCyto100300.sequence{24};
seqB = yCyto100300.sequence{47};
numStaples = 10;
numTrials = 500;
results.AA = staplingsEnergyDistribution1Pair(seqA,seqA,numStaples,numTrials);
results.AB = staplingsEnergyDistribution1Pair(seqA,seqB,numStaples,numTrials);
end