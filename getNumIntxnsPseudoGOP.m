function [results] = getNumIntxnsPseudoGOP()

load phycfull;
load godata;
%load y2hdata;

results.numProt = nan(100,1);
results.allStats = nan(100,3);
results.selfOn = nan(100,1);
%results.sameStats = zeros(100,3);

for i=1:length(godata.pIdx)
    allPossPairs = nchoosek(godata.pIdx{i},2);
    if length(allPossPairs) == 1
        continue
    end
   [iOn,iOff,iTotal,selfOn,iOnSame,iOffSame,iSameAll] = countNumOnOffIntrxns(allPossPairs,godata.pIdx{i},phycfull); 
    results.allStats(i,:) = [iOn iOff iTotal];
    results.selfOn(i) = length(selfOn);
    %[r,~] = size(iSameAll);
    %results.sameStats(i,:) = [iOnSame,iOffSame,r];
    results.numProt(i) = length(godata.pIdx{i});
end


end