function [results] = runIsingHomologs(nTrials)

% Load organism data
%addpath('C:\Users\Kelly\Dropbox (MIT)\EnglandRotation\phCrossSpecies');
addpath('phCrossSpecies');
results.speciesNames = {'acarolinensis','afumigatus','agossypii','amexicanus',...
    'aniger','aoryzae','btaurus','celegans2','cintestinalis','cneoformans'...
    'drerio','drosophila','ggallus','human','kpastoris','lchalumnae',...
    'mmusculus','pmarinus','spombe','tguttata','tmelanosporum','ttruncatus',...
    'ylipolytica','yeast'};

N = length(results.speciesNames);
results.data = cell(N,1);
results.isingData = cell(N,1);
results.percentLowerDE = zeros(N,1);
results.percentLowerAcidic = zeros(N,1);
results.percentLowerNeutral = zeros(N,1);
results.numIntrxns = zeros(N,1);
results.actDE = zeros(N,1);
results.actNeutral = zeros(N,1);
results.actAcidic = zeros(N,1);

results.isingDataOS = cell(N,1);
results.percentLowerDEOS = zeros(N,1);
results.percentLowerAcidicOS = zeros(N,1);
results.percentLowerNeutralOS = zeros(N,1);

results.actDEOS = zeros(N,1);
results.actNeutralOS = zeros(N,1);
results.actAcidicOS = zeros(N,1);

for i=1:N
    disp(i);
    results.data{i} = importdata(strcat(results.speciesNames{i},'Y2H.mat'));
    
    % Add in sequence length
    results.data{i}.seqLength = zeros(length(results.data{i}.seq),1);
    for j=1:length(results.data{i}.seq)
       
        results.data{i}.seqLength(j) = length(results.data{i}.seq{j});
        if results.data{i}.seq{j}(end) == '*'
            results.data{i}.seqLength(j) = results.data{i}.seqLength(j)-1;
            results.data{i}.seq{j}(end) = [];
        end
    end
    
    % Get results using purely yeast sequences
    results.isingData{i} = getIsingHist(nTrials,results.data{i}.intrxnsBasic);
    results.percentLowerAcidic(i) = length(find(results.isingData{i}.eTotalAcidic < results.isingData{i}.actETotalAcidic))/nTrials;
    results.percentLowerNeutral(i) = length(find(results.isingData{i}.eTotalNeutral < results.isingData{i}.actETotalNeutral))/nTrials;
    results.percentLowerDE(i) = length(find(results.isingData{i}.dEHist < results.isingData{i}.actdE))/nTrials;
    [results.numIntrxns(i),~] = size(results.data{i}.intrxnsBasic);
    results.actDE(i) = results.isingData{i}.actdE;
    results.actNeutral(i) = results.isingData{i}.actETotalNeutral;
    results.actAcidic(i) = results.isingData{i}.actETotalAcidic;
    
    % Get results using organism (cross-species) sequences
    results.isingDataOS{i} = getIsingHist(nTrials,results.data{i}.intrxnsBasHID, results.data{i});
    
    results.percentLowerAcidicOS(i) = length(find(results.isingDataOS{i}.eTotalAcidic < results.isingDataOS{i}.actETotalAcidic))/nTrials;
    results.percentLowerNeutralOS(i) = length(find(results.isingDataOS{i}.eTotalNeutral < results.isingDataOS{i}.actETotalNeutral))/nTrials;
    results.percentLowerDEOS(i) = length(find(results.isingDataOS{i}.dEHist < results.isingDataOS{i}.actdE))/nTrials;
    %[results.numIntrxnsOS(i),~] = size(results.data{i}.intrxnsBasic);
    results.actDEOS(i) = results.isingDataOS{i}.actdE;
    results.actNeutralOS(i) = results.isingDataOS{i}.actETotalNeutral;
    results.actAcidicOS(i) = results.isingDataOS{i}.actETotalAcidic;

end

