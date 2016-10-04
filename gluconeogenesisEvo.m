function [results] = gluconeogenesisEvo()

addpath('phCrossSpecies');
prots = [163 428];

results.speciesNames = {'acarolinensis','afumigatus','agossypii','amexicanus',...
    'aniger','aoryzae','btaurus','celegans2','cintestinalis','cneoformans'...
    'drerio','drosophila','ggallus','human','kpastoris','lchalumnae',...
    'mmusculus','pmarinus','spombe','tguttata','tmelanosporum','ttruncatus',...
    'ylipolytica','yeast'};
nSpecies = length(results.speciesNames);

results.data = cell(nSpecies,1);
results.eNeutralAB = nan(nSpecies,1);
results.eAcidicAB = nan(nSpecies,1);
results.dEAB = nan(nSpecies,1);

results.eNeutralAA = nan(nSpecies,1);
results.eAcidicAA = nan(nSpecies,1);
results.dEAA = nan(nSpecies,1);

results.eNeutralBB = nan(nSpecies,1);
results.eAcidicBB = nan(nSpecies,1);
results.dEBB = nan(nSpecies,1);

for i=1:nSpecies
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
    
    hprots = find(ismember(results.data{i}.yIdx,prots));
    if length(hprots) == 2
        hprotsPairs = [hprots(1) hprots(2); hprots(1) hprots(1); hprots(2) hprots(2)];
        hprotsPairs
        r = isingModelIntrxns(hprotsPairs,7,5,results.data{i});
        
        
        results.eNeutralAB(i) = r.eComponentsNeutral(1);
        results.eAcidicAB(i) = r.eComponentsAcidic(1);
        results.dEAB(i) = r.eComponentsAcidic(1)-r.eComponentsNeutral(1);
        
        results.eNeutralAA(i) = r.eComponentsNeutral(2);
        results.eAcidicAA(i) = r.eComponentsAcidic(2);
        results.dEAA(i) = r.eComponentsAcidic(2)-r.eComponentsNeutral(2);
        
        results.eNeutralBB(i) = r.eComponentsNeutral(3);
        results.eAcidicBB(i) = r.eComponentsAcidic(3);
        results.dEBB(i) = r.eComponentsAcidic(3)-r.eComponentsNeutral(3);
    end
    
end

nonspecies = find(isnan(results.eNeutralAB));
results.eNeutralAB(nonspecies) = [];
results.eAcidicAB(nonspecies) = [];
results.dEAB(nonspecies) = [];

results.eNeutralAA(nonspecies) = [];
results.eAcidicAA(nonspecies) = [];
results.dEAA(nonspecies) = [];

results.eNeutralBB(nonspecies) = [];
results.eAcidicBB(nonspecies) = [];
results.dEBB(nonspecies) = [];
results.speciesNames(nonspecies) = [];

end