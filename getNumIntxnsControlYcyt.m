function [results] =  getNumIntxnsControlYcyt(nIts,nRTS)

% nIts is number of iterations to go through
load phycfull;
load y2hdata;
load yCyt;

% Get random assortments of 28 proteins in yCyt
%nRTS = 28;

pIdx = zeros(nIts,nRTS);
z = 1:length(yCyt.id);

for i=1:nIts
    z2 = randperm(length(z));
    pIdx(i,:) = z2(1:nRTS);
end

results.numProt = nRTS * ones(nIts,1);
results.allStats = zeros(nIts,3);
results.selfOn = zeros(nIts,1);
results.sameStats = zeros(nIts,3);

for i=1:nIts
    
   [iOn,iOff,iTotal,selfOn,iOnSame,iOffSame,iSameAll] = countNumOnOffIntrxns(y2hdata.idxCytPairs,pIdx(i,:)',phycfull); 
    results.allStats(i,:) = [iOn iOff iTotal];
    results.selfOn(i) = length(selfOn);
    [r,~] = size(iSameAll);
    results.sameStats(i,:) = [iOnSame,iOffSame,r];
    %results.numProt(i) = length(godata.pIdx{i});
end
end