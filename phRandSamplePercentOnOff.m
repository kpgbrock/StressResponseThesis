function [results] = phRandSamplePercentOnOff(nTrials)

load y2hdata;
load phycfull;

ly2h = length(unique([y2hdata.idxNonSelf(:,1);y2hdata.idxNonSelf(:,2)]));
ny2hInt = length(y2hdata.idxNonSelf);

results.numOnToOff = zeros(nTrials,1);
results.numStayOn = zeros(nTrials,1);
results.numOn7 = zeros(nTrials,1);
results.numOffToOn = zeros(nTrials,1);
results.numStayOff = zeros(nTrials,1);
results.numOff7 = zeros(nTrials,1);

for i=1:nTrials
    lProts = randperm(1804);
    lProts = lProts(1:ly2h);
    z = nchoosek(lProts,2);
    x = randperm(length(z));
    x = x(1:ny2hInt);
    
    lInt = z(x,:);
    [temp] = phPercentageOnPredicated(lInt);
    results.numOnToOff(i) = temp.numOnToOff;
	results.numStayOn(i) = temp.numStayOn;
    results.numOn7(i) = temp.numOn7;
    results.numOffToOn(i) = temp.numOffToOn;
    results.numStayOff(i) = temp.numStayOff;
    results.numOff7(i) = temp.numOff7;
end

% results.percentOnToOff = results.numOnToOff./results.numOn7;
% results.percentOffToOn(i) = results.numOffToOn./results.numOff7;
end