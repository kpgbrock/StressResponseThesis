function [results] = glycProtAnalysis(nTrials)

load glycProt;
load yCyt;
load phycfull;

N = length(glycProt.locsTrunc);
protPairs = nchoosek(glycProt.locsTrunc,2);

resultsAct = isingModelIntrxns(protPairs, 7, 5);
results.actETotalNeutral = resultsAct.eTotalNeutral;
results.actETotalAcidic = resultsAct.eTotalAcidic;
results.actDE = resultsAct.dE;

results.eTotalNeutral = zeros(nTrials,1);
results.eTotalAcidic = zeros(nTrials,1);
results.dE = zeros(nTrials,1);

for i=1:nTrials
    %disp(i);
    z = randperm(length(yCyt.id),N);
    protPairsRand = nchoosek(z,2);
    
    z2 = isingModelIntrxns(protPairsRand, 7, 5);
    results.eTotalNeutral(i) = z2.eTotalNeutral;
    results.eTotalAcidic(i) = z2.eTotalAcidic;
    results.dE(i) = z2.dE;
end
    

end