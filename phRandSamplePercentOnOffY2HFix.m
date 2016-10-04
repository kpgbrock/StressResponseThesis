function [results] = phRandSamplePercentOnOffY2HFix(nTrials)

load y2hdata;
load phycfull;

% setOfIntrxns = nchoosek(y2hdata.idxIndProtAll,2);
setOfIntrxns = nchoosek(y2hdata.idxIndProtNS,2);
numPossInt = length(setOfIntrxns);

results.numOnToOff = zeros(nTrials,1);
results.numStayOn = zeros(nTrials,1);
results.numOn7 = zeros(nTrials,1);
results.numOffToOn = zeros(nTrials,1);
results.numStayOff = zeros(nTrials,1);
results.numOff7 = zeros(nTrials,1);

for i=1:nTrials
    z = randperm(numPossInt);
    %z = z(1:length(y2hdata.idxCytPairs));
    z = z(1:length(y2hdata.idxNonSelf));
    if mod(i,500) == 0
        disp(i);
    end
    lInt = setOfIntrxns(z,:);
    [temp] = phPercentageOnPredicatedGeneral(lInt,phycfull);
    results.numOnToOff(i) = temp.numOnToOff;
	results.numStayOn(i) = temp.numStayOn;
    results.numOn7(i) = temp.numOn7;
    results.numOffToOn(i) = temp.numOffToOn;
    results.numStayOff(i) = temp.numStayOff;
    results.numOff7(i) = temp.numOff7;
end

results.percentOnToOff = results.numOnToOff./results.numOn7;
results.percentOffToOn = results.numOffToOn./results.numOff7;
end