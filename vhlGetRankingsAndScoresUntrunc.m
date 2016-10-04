load upVHLPairSwapAll;
v = upVHLPairSwapAll;
load vhlUP;
ps = vhlUP.pairSwaps;

r = rankOrderPairSwapVHLUntrunc;

results = zeros(20,2);
for i=1:20
    results(i,1) = find((r.sumScoreRanked(:,3) == min(ps(i,:))) & (r.sumScoreRanked(:,4) == max(ps(i,:))));
    results(i,2) = r.sumScoreRanked(results(i,1),1);
end