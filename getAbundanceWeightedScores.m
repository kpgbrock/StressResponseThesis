function [results] = getAbundanceWeightedScores(nTrials)

load phyce;
load yCytoExp;

load y2hdata;
load phycfull;


N = length(yCytoExp.abundance);
results.phyceProt = zeros(nTrials,2);
results.phyceScore5 = zeros(nTrials,1);
results.phyceScore7 = zeros(nTrials,1);

z5 = find(phyce.pHRange == 5);
z7 = find(phyce.pHRange == 7);

for i=1:nTrials
    p1 = randsample(N,1,'true',yCytoExp.abundance/sum(yCytoExp.abundance));
    p2 = randsample(N,1,'true',yCytoExp.abundance/sum(yCytoExp.abundance));
    
    results.phyceProt(i,:) = [p1 p2];
    results.phyceScore5(i) = phyce.charge(z5,p1)*phyce.charge(z5,p2)/(phyce.seqLength(p1)*phyce.seqLength(p2));
    results.phyceScore7(i) = phyce.charge(z7,p1)*phyce.charge(z7,p2)/(phyce.seqLength(p1)*phyce.seqLength(p2));
end

end
