function [results] = phPercentageOnPredicatedGeneral(intrxns,data)

%load y2hdata;
%load phycfull;

z5 = find(data.pHRange == 5);
z7 = find(data.pHRange == 7);

% results.findAlreadyOn5 = 0;
% results.findAlreadyOff5 = 0;
findAlreadyOn7 = 0;
findAlreadyOff7 = 0;

z = [intrxns data.charge(z5,intrxns(:,1))' data.charge(z5,intrxns(:,2))' data.charge(z7,intrxns(:,1))' data.charge(z7,intrxns(:,2))'];

% results.findAlreadyOn5 = find(z(:,3).*z(:,4) < 0);
% results.findAlreadyOff5 = find(z(:,3).*z(:,4) > 0);
findAlreadyOn7 = find(z(:,5).*z(:,6) < 0);
findAlreadyOff7 = find(z(:,5).*z(:,6) > 0);

[testOn,iOff] = getOnOffIntrxns(intrxns(findAlreadyOn7,:),data);

[iOn,testOff] = getOnOffIntrxns(intrxns(findAlreadyOff7,:),data);
disp(testOn);
disp(testOff);

[r,c] = size(iOff);
results.numOnToOff = r;
results.numStayOn = length(findAlreadyOn7) - results.numOnToOff;

[r,c] = size(iOn);
results.numOn7 = length(findAlreadyOn7);
results.numOffToOn = r;
results.numStayOff = length(findAlreadyOff7) - results.numOffToOn;
results.numOff7 = length(findAlreadyOff7);
end