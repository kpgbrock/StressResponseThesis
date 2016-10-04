function [results] = getGOWithPuncta()

load yCyt;
load phycfull;
load godata;
load ellingtonS4;
load puncta180;

results.statsP180 = zeros(100,3);
results.statsEll = zeros(100,3);
results.p180NumSame = zeros(100,1);
results.ellNumSame = zeros(100,1);

ngo = length(godata.pIdx);

% Get charges for interacting pairs at 5 and 7
z5 = find(phycfull.pHRange == 5);
z7 = find(phycfull.pHRange == 7);

for i=1:ngo
   prot1 = godata.pIdx{i};
   prot2 = puncta180.idxCytFull;
   prot3 = ellingtonS4.pIdx;

   [pairsPuncta,xP] = getPossIntTwoGroups(prot1,prot2);
   [pairsCyc,xC] = getPossIntTwoGroups(prot1,prot3);
   
   
   %
   results.p180NumSame(i) = xP;
   results.ellNumSame(i) = xC;  
   
   pairs = pairsPuncta;
   [r,~] = size(pairs);
   z = [phycfull.charge(z5,pairs(:,1))' phycfull.charge(z5,pairs(:,2))' phycfull.charge(z7,pairs(:,1))' phycfull.charge(z7,pairs(:,2))'];
   zOn = pairs((z(:,1).*z(:,2) < 0) & (z(:,3).*z(:,4) > 0),:);
   zOff = pairs((z(:,1).*z(:,2) > 0) & (z(:,3).*z(:,4) < 0),:);
   [r1,~] = size(zOn);
   [r2,~] = size(zOff);
   results.statsP180(i,:) = [r1 r2 r];
   
   pairs = pairsCyc;
   [r,~] = size(pairs);
   z = [phycfull.charge(z5,pairs(:,1))' phycfull.charge(z5,pairs(:,2))' phycfull.charge(z7,pairs(:,1))' phycfull.charge(z7,pairs(:,2))'];
   zOn = pairs((z(:,1).*z(:,2) < 0) & (z(:,3).*z(:,4) > 0),:);
   zOff = pairs((z(:,1).*z(:,2) > 0) & (z(:,3).*z(:,4) < 0),:);
   [r1,~] = size(zOn);
   [r2,~] = size(zOff);
   results.statsEll(i,:) = [r1 r2 r];
end

results.statsP180(:,4) = results.statsP180(:,1)./results.statsP180(:,3);
results.statsEll(:,4) = results.statsEll(:,1)./results.statsEll(:,3);

end