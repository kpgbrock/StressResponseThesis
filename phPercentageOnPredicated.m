function [results] = phPercentageOnPredicated(intrxns)

%load y2hdata;
load phycfull;

z5 = find(phycfull.pHRange == 5);
z7 = find(phycfull.pHRange == 7);

% results.findAlreadyOn5 = 0;
% results.findAlreadyOff5 = 0;
findAlreadyOn7 = 0;
findAlreadyOff7 = 0;

z = [intrxns phycfull.charge(z5,intrxns(:,1))' phycfull.charge(z5,intrxns(:,2))' phycfull.charge(z7,intrxns(:,1))' phycfull.charge(z7,intrxns(:,2))'];

% results.findAlreadyOn5 = find(z(:,3).*z(:,4) < 0);
% results.findAlreadyOff5 = find(z(:,3).*z(:,4) > 0);
findAlreadyOn7 = find(z(:,5).*z(:,6) < 0);
findAlreadyOff7 = find(z(:,5).*z(:,6) > 0);

[testOn,iOff] = getOnOffIntrxns(intrxns(findAlreadyOn7,:),phycfull);
[iOn,testOff] = getOnOffIntrxns(intrxns(findAlreadyOff7,:),phycfull);
disp(testOn);
disp(testOff);
results.numOnToOff = length(iOff);
results.numStayOn = length(findAlreadyOn7) - results.numOnToOff;
results.numOn7 = length(findAlreadyOn7);
results.numOffToOn = length(iOn);
results.numStayOff = length(findAlreadyOff7) - results.numOffToOn;
results.numOff7 = length(findAlreadyOff7);
end