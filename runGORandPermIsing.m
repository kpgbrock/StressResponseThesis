function [results] = runGORandPermIsing(nTrials)

load yCyt;
load godata;
load phycfull;

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
        
        seqString = strjoin(yCyt.sequence(godata.pIdx{i}));
        seqString(seqString == '*') = [];
        seqString(seqString == ' ') = [];
        seqLengths = phycfull.seqLength(godata.pIdx{i});
        %disp(seqLengths);
        % Compare to random drawings
        for j=1:nTrials
            seqString = seqString(randperm(length(seqString)));
            tempdata.seq = mat2cell(seqString,1,seqLengths);
            tempdata.seqLength = seqLengths;
            tempdata.pHRange = phycfull.pHRange;
            tempdata.charge = zeros(length(phycfull.pHRange),length(tempdata.seq));
            
            % Get charge array:  findPICurveWiki
            for k=1:length(tempdata.seq)
                [fcharge,~] = findPICurveWiki(tempdata.seq{k},tempdata.pHRange);
                tempdata.charge(:,k) = fcharge;
            end
            %
            
            randIntrxns = nchoosek(1:length(tempdata.seq),2);
            
            
            data2 = isingModelIntrxns(randIntrxns,7,5, tempdata);
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