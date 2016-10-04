function [results] = vhlGetRankingsAndScoresAllMuts()
load vhlUP;
load upVHLPointMutResults;
load upVHLPairSwapResults;

r = rankOrderPairSwapVHL;
r2 = rankOrderPointMutVHL;

x = 63:204;

results = nan(20,3);
for i=1:20
    p1 = vhlUP.pairSwaps(i,1);
    p2 = vhlUP.pairSwaps(i,2);
    i
    
    if (any(x==p1) && any(x==p2))
        results(i,1) = find((r.sumScoreRanked(:,3) == min([p1 p2])) & (r.sumScoreRanked(:,4) == max([p1 p2])));
        results(i,2) = r.sumScoreRanked(results(i,1),1);
        results(i,3) = 2;
    elseif (any(x==p1))
        results(i,1) = find((r2.sumScoreRanked(:,3) == p1) & (r2.sumScoreRanked(:,4) == aa2int(vhlUP.orgSeq(p2))));
        results(i,2) = r2.sumScoreRanked(results(i,1),1);
        results(i,3) = 1;
    elseif any(x==p2)
        results(i,1) = find((r2.sumScoreRanked(:,3) == p2) & (r2.sumScoreRanked(:,4) == aa2int(vhlUP.orgSeq(p1))));
        results(i,2) = r2.sumScoreRanked(results(i,1),1);
        results(i,3) = 1;
    end
end
end