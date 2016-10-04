function [results] = getNumIntxnsAllGOP()

load phycfull;
load godata;
load y2hdata;

results.numProt = zeros(100,1);
results.allStats = zeros(100,3);
results.selfOn = zeros(100,1);
results.sameStats = zeros(100,3);

for i=1:length(godata.pIdx)
    
   [iOn,iOff,iTotal,selfOn,iOnSame,iOffSame,iSameAll] = countNumOnOffIntrxns(y2hdata.idxCytPairs,godata.pIdx{i},phycfull); 
    results.allStats(i,:) = [iOn iOff iTotal];
    results.selfOn(i) = length(selfOn);
    [r,~] = size(iSameAll);
    results.sameStats(i,:) = [iOnSame,iOffSame,r];
    results.numProt(i) = length(godata.pIdx{i});
end
end