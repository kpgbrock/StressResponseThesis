function [results] = phRandSampleDiffChProfileGen(lProts,nTrials,phLow,phHi)

load phycfull;

N = 1804;
npH = length(phycfull.pHRange);

results.phL = find(phycfull.pHRange == phLow);
results.phH = find(phycfull.pHRange == phHi);

results.avgSelfIntCh = mean(phycfull.charge(:,lProts)')';
results.avgAllCh = mean(phycfull.charge(:,:)')';
results.avgSelfIntLength = mean(phycfull.chPerLen(:,lProts)')';
results.avgAllChLength = mean(phycfull.chPerLen(:,:)')';

% What we want to compare our nTrials values to
results.magSelfIntCh = sum(abs(results.avgSelfIntCh(results.phL:results.phH) - results.avgAllCh(results.phL:results.phH)));
results.magSelfIntChLen = sum(abs(results.avgSelfIntLength(results.phL:results.phH) - results.avgAllChLength(results.phL:results.phH)));

results.avgChTrials = zeros(npH,nTrials);
results.avgChLenTrials = zeros(npH,nTrials);
results.stdChTrials = zeros(npH,nTrials);
results.stdChTrials = zeros(npH,nTrials);
results.magDiffCh = zeros(nTrials,1);
results.magDiffChLen = zeros(nTrials,1);


for i=1:nTrials
   z = randperm(N);
   z = z(1:length(lProts));
   
   newCh = phycfull.charge(:,z);
   newChLen = phycfull.chPerLen(:,z);
%    size(mean(newCh')')
%    size(results.avgChTrials(:,i))
   results.avgChTrials(:,i) = mean(newCh')';
   results.avgChLenTrials(:,i) = mean(newChLen')';
   results.stdChTrials(:,i) = std(newCh')';
   results.stdChLenTrials(:,i) = std(newChLen')';
   
   results.magDiffCh(i) = sum(abs(results.avgChTrials(results.phL:results.phH,i) - results.avgAllCh(results.phL:results.phH)));
   results.magDiffChLen(i) = sum(abs(results.avgChLenTrials(results.phL:results.phH,i) - results.avgAllChLength(results.phL:results.phH)));
   
   
end
end