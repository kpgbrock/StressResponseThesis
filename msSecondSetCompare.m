function [results] = msSecondSetCompare(nTrials)

load yCyt;
load phycfull;
load yMassSpec;

% Get list of real proteins that interact with tsr2
tsr2ID = 1227;
x7 = [yMassSpec.id7; tsr2ID];
N = length(x7);

% Find where ph 5 and 7 are located in data structure
ph5 = find(phycfull.pHRange == 5);
ph7 = find(phycfull.pHRange == 7);

% Calculate q1*q2/(L1*L2) for protein pairings
[results.ch5, results.ch7, results.actSum5, results.actSum7, results.actDE, results.actSumDE] = getChargePairings(x7,ph5, ph7, phycfull);
[results.tsrCh5, results.tsrCh7, results.tsrActSum5, results.tsrActSum7, results.tsrActDE, results.tsrActSumDE] = getChargePairings([x7(1:end-1) tsr2ID*ones(length(x7)-1,1)],ph5,ph7,phycfull);
results.distSum5 = zeros(nTrials,1);
results.distSum7 = zeros(nTrials,1);
results.distDE = zeros(nTrials,1);

for i=1:nTrials
   
   % Get random protein ids that don't include tsr2ID
   xNew = [1:(tsr2ID-1) tsr2ID+1:length(yCyt.id)]; 
   
   xNew = xNew(randperm(length(xNew),N-1))';
   
   [~,~, results.tsrDistSum5(i), results.tsrDistSum7(i), ~, results.tsrDistDE(i)] = getChargePairings([xNew tsr2ID*ones(length(xNew),1)],ph5, ph7, phycfull);
   
   xNew = [xNew; tsr2ID];
   
   [~,~, results.distSum5(i), results.distSum7(i), ~, results.distDE(i)] = getChargePairings(xNew,ph5, ph7, phycfull);

end
   

end

function [dist5, dist7, sum5, sum7, dE, sumDE] = getChargePairings(protList,ph5, ph7, phycfull)

[~,c] = size(protList);
if c == 2
    z = protList;
else
    z = nchoosek(protList,2);
end

dist5 = phycfull.chPerLen(ph5,z(:,1)).*phycfull.chPerLen(ph5,z(:,2));
dist7 = phycfull.chPerLen(ph7,z(:,1)).*phycfull.chPerLen(ph7,z(:,2));
sum5 = sum(dist5);
sum7 = sum(dist7);
dE = dist5 - dist7;
sumDE = sum(dE);

end