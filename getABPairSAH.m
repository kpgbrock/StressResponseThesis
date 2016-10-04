function [results] = getABPairSAH()
addpath('..');
load yCyto50;
yc = yCyto50;

n = 50;
l = (n*(n-1)/2);
results.AB = zeros(l,3);
results.AA = zeros(n,2);
ct = 1;
for i=1:50
    disp(i);
    for j=1:(i-1)
        [~,~,~,emin,~,~,~,~] = pholderAggSAH(yc.sequence{i},yc.sequence{j});
        results.AB(ct,1) = i;
        results.AB(ct,2) = j;
        results.AB(ct,3) = emin;
        ct = ct+1;
    end
    [~,~,~,emin,~,~,~,~] = pholderAggSAH(yc.sequence{i},yc.sequence{i});
    results.AA(i,1) = i;
    results.AA(i,2) = emin;
end

end

