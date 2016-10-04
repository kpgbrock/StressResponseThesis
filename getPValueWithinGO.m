function [results] = getPValueWithinGO(nIts)

load godata;

% Get how many proteins are involved in each GO group
N = length(godata.pIdx);
nL = zeros(N,1);
for i=1:N
   nL(i) = length(godata.pIdx{i});
end
nP = unique(nL);

N2 = length(nP);
results.nP = nP;
results.dist = cell(N2,1);
results.percentile = nan(N,1);

% For each possible GO group size, do nIts random trials involving x = size
% of GO group proteins
for i=2:N2
   disp(i);
   tempres = getNumIntxnsControlYcytAllPossPairs(nIts,nP(i));
   tempres.allStats(:,4) = tempres.allStats(:,1)./tempres.allStats(:,3);
   results.dist{i} = tempres.allStats(:,4);   
end

% Get fraction of on interactions for actual GO groups
[actRes] = getNumIntxnsPseudoGOP();
actRes.allStats(:,4) = actRes.allStats(:,1)./actRes.allStats(:,3);

% Calculate percentiles based on our trials and our actual results
for i=1:N    
    if nL(i)==1
        continue;
    end
    d = results.dist{nP == nL(i)};     
    results.percentile(i) = length(find(d < actRes.allStats(i,4)))/nIts;
end
end