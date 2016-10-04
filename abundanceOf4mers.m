function [results] = abundanceOf4mers()

load yCytoExp;
yc = yCytoExp;

%yc.abundance = yc.abundance(randperm(length(yc.abundance)));

results = zeros(20,20,20,20);

for i=1:length(yc.abundance)
   seq = yc.sequence{i};
   seq(end) = [];
   seq = seq(randperm(length(seq)));
   domains = breakIntoDomains(seq,4);
   domains(end) = [];
   
   for j=1:length(domains)
      domainInt = aa2int(domains{j});
      temp = results(domainInt(1),domainInt(2),domainInt(3),domainInt(4));
      results(domainInt(1),domainInt(2),domainInt(3),domainInt(4)) = temp + yc.abundance(i);
   end
end
end