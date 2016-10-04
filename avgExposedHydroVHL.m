function [results] = avgExposedHydroVHL(cutoff)

load vhlUP;
xrange = 64:203;
results = zeros(21,4);
results(:,1) = (0:20)';

for i=1:21
    bt = vhlUP.btTrunc(:,i);
    bt = bt/max(bt);
    
    if i == 1
        seq = vhlUP.orgSeq(xrange);
    else
        seq = vhlUP.mutSeq{i-1}(xrange);
    end
    
    n = find(bt(bt > cutoff));
    hydro = zeros(length(n),1);
    for j=1:n
        hydro(j) = KD(seq(n(j)));
    end
    results(i,2) = mean(hydro);
    results(i,3) = std(hydro);
    results(i,4) = length(n);
end