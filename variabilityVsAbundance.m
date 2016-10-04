function [results] = variabilityVsAbundance

load yCytoExp;
results = zeros(456,2);
ct = 1;
for i=1:length(yCytoExp.id)
    
   seq = yCytoExp.sequence{i};
   if (length(seq) >= 100) && (length(seq) <= 300)
       
       [score,~,~] = pholderVariance(seq,25);
       results(ct,:) = [yCytoExp.abundance(i),score];
       ct = ct+1;
       disp(i);
   end
end
plot(results(:,1),results(:,2),'*');
end