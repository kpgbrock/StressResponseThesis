function [results] = goPiePlotMHMC()

load yCyt;
load godata;
load yeastGOAnnotation;

%listGO = {'generation of precursor metabolites and energy'; 'mitotic cell cycle'; 'meiotic cell cycle'; 'ribosomal small subunit biogenesis'; 'regulation of cell cycle'; 'cytoplasmic translation'; 'RNA modification'; 'carbohydrate metabolic process'; 'other'};
listGO = {'mitotic cell cycle'; 'meiotic cell cycle'; 'regulation of cell cycle'; 'ribosomal small subunit biogenesis'; 'cytoplasmic translation'; 'RNA modification'; 'carbohydrate metabolic process'; 'generation of precursor metabolites and energy'; 'other'};
intProt = [51;153;208;211;221;228;243;290;295;317;357;361;428;472;548;557;570;580;644;820;903;1017;1025;1227;1272;1344;1379;1386;1518;1745];
protAll = 1:1804;

[pvalsMC,multiCtMC] = goPiePlotHelper(intProt,listGO(1:end-1),godata,yeastGOAnnotation);
[pvalsAll,multiCtAll] = goPiePlotHelper(protAll,listGO(1:end-1),godata,yeastGOAnnotation);

results.pvalsMC = pvalsMC;
results.pvalsAll = pvalsAll;
results.multiCtMC = multiCtMC;
results.multiCtAll = multiCtAll;

subplot(1,2,1);
%pie(pvalsMC,listGO);
pie(pvalsMC);
%title('Monte Carlo N = 30');
subplot(1,2,2);
%pie(pvalsAll,listGO);
pie(pvalsAll);
%title('All 1804 Cytosolic Proteins');

end

function [pvals,multCt] = goPiePlotHelper(protList,listGO,godata,yeastGOAnnotation)

pvals = zeros(length(listGO),1);
multCt = zeros(length(protList),1);
for i=1:length(protList)
    
    for j=1:length(listGO)
        %yeastGOAnnotation.anno(godata.goIdxP{j})
        %listGO{i}
        if ~isempty(find(ismember(yeastGOAnnotation.anno(godata.goIdxP{protList(i)}),listGO{j}), 1))
            pvals(j) = pvals(j) + 1;
            multCt(i) = multCt(i) + 1;
        end
    end
    
end

pvals = [pvals;length(multCt(multCt == 0))];
end

