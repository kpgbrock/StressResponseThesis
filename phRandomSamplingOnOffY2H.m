function [results] = phRandomSamplingOnOffY2H(nTrials)

load y2hdata;
load phycfull;

N = length(y2hdata.idxNonSelf);

z5 = find(phycfull.pHRange == 5);
z7 = find(phycfull.pHRange == 7);

p = nchoosek((1:1804)',2);
x = length(p);
results.iOn = zeros(nTrials,1);
results.iOff = zeros(nTrials,1);

for i=1:nTrials
   z = randperm(x);
   z=z(1:N);
   
   
   [iOn,iOff] = getOnOffIntrxns(p(z,:),phycfull);
   [r,~] = size(iOn);
   results.iOn(i) = r/N;
   [r,~] = size(iOff);
   results.iOff(i) = r/N;
end

end