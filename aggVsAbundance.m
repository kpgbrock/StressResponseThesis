function [results] = aggVsAbundance()

load yCytoExp;
results = nan(length(yCytoExp.abundance),2);
for i=1:length(yCytoExp.abundance)
    seq = yCytoExp.sequence{i};
    if (length(seq) >= 100) && (length(seq) <= 300)
       [~,~,~,emin,~,~,~,~] = pholderAgg(seq,seq);
       results(i,1) = yCytoExp.abundance(i);
       results(i,2) = emin;
    end
end
end