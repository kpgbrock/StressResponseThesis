function [results] = runGOGenericIsing(nTrials)

load yCyt;
load godata;

N = length(godata.pIdx);
results.numProts = zeros(N,1);
results.numIntrxns = zeros(N,1);
results.actAcidic = nan(N,1);
results.actNeutral = nan(N,1);
results.actDE = nan(N,1);
results.percentLowerAcidic = nan(N,1);
results.percentLowerNeutral = nan(N,1);
results.percentLowerDE = nan(N,1);
results.histAcidic = nan(nTrials,N);
results.histNeutral = nan(nTrials,N);
results.histDE = nan(nTrials,N);

for i=1:N
    
    disp(i);
    % Get energy scores for all possible pairwise interactions within GO
    % groups
    
    results.numProts(i) = length(godata.pIdx{i});
    if results.numProts(i) >= 2
        protIntrxns = nchoosek(godata.pIdx{i},2);
        [results.numIntrxns(i),~] = size(protIntrxns);
        
        [data] = isingModelIntrxns(protIntrxns, 7,5);
        results.actAcidic(i) = data.eTotalAcidic;
        results.actNeutral(i) = data.eTotalNeutral;
        results.actDE(i) = data.dE;
        
        % Compare to random drawings
        for j=1:nTrials
            randProt = randperm(length(yCyt.id),results.numProts(i));
            randIntrxns = nchoosek(randProt,2);
            data2 = isingModelIntrxns(randIntrxns,7,5);
            results.histAcidic(j,i) = data2.eTotalAcidic;
            results.histNeutral(j,i) = data2.eTotalNeutral;
            results.histDE(j,i) = data2.dE;
        end
        
        results.percentLowerAcidic(i) = length(find(results.histAcidic(:,i) < results.actAcidic(i)))/nTrials;
        results.percentLowerNeutral(i) = length(find(results.histNeutral(:,i) < results.actNeutral(i)))/nTrials;
        results.percentLowerDE(i) = length(find(results.histDE(:,i) < results.actDE(i)))/nTrials;
    end
    
end
end