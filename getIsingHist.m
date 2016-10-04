function [results] = getIsingHist(nTrials,protPairs,phycfull)

%load phycfull;
if nargin==1
    %load y2hdata;    
    protPairs = y2hdata.idxNonSelf;
elseif nargin == 2
    load phycfull;
end

% New line
protIDs = unique([protPairs(:,1); protPairs(:,2)]);
%protIDs = y2hdata.idxIndProtNS;
%[numPairs,~] = size(y2hdata.idxNonSelf);
[numPairs,~] = size(protPairs);

idPairs = nchoosek(protIDs,2);

[nPoss,~] = size(idPairs);

results.dEHist = zeros(nTrials,1);
%results.eComponentsNeutral = zeros(nTrials,1);
results.eTotalNeutral = zeros(nTrials,1);
results.eTotalAcidic = zeros(nTrials,1);
results.eCompNeutral = cell(nTrials,1);
results.eCompAcidic = cell(nTrials,1);

%actResults = isingModelIntrxns(y2hdata.idxNonSelf,5,7);
%actResults = isingModelIntrxns(y2hdata.idxNonSelf,7,5);

% actResults = isingModelIntrxns(protPairs,7,5);
actResults = isingModelIntrxns(protPairs,7,5,phycfull);

results.actdE = actResults.dE;
results.actETotalNeutral = actResults.eTotalNeutral;
results.actETotalAcidic = actResults.eTotalAcidic;
results.actEComponentsNeutral = actResults.eComponentsNeutral;
results.actEComponentsAcidic = actResults.eComponentsAcidic;

% To store interactome trials
results.trials = cell(nTrials,1);


for i=1:nTrials
   
    trialPairs = idPairs(datasample(1:nPoss,numPairs,'Replace',false),:);
    results.trials{i} = trialPairs;
    
    %tempResults = isingModelIntrxns(trialPairs,7,5);
    tempResults = isingModelIntrxns(trialPairs,7,5,phycfull);
    
    results.dEHist(i) = tempResults.dE;
    results.eTotalNeutral(i) = tempResults.eTotalNeutral;
    results.eTotalAcidic(i) = tempResults.eTotalAcidic;
    results.eCompNeutral{i} = tempResults.eComponentsNeutral;
    results.eCompAcidic{i} = tempResults.eComponentsAcidic;
end
end