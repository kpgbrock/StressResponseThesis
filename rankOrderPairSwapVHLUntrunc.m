function [results] = rankOrderPairSwapVHLUntrunc()

load upVHLPairSwapAll;
upVHLPairSwapResults = upVHLPairSwapAll;

% Normalize the burial traces
bt = upVHLPairSwapResults.bt;
N = length(upVHLPairSwapResults.pairSwap);
for i=1:N
    bt(:,i) = bt(:,i)/max(bt(:,i));
end

% Get actual pair swap positions before truncating
ps = upVHLPairSwapResults.pairSwap;

results.box1min = [ps zeros(N,1)];
results.box2min = [ps zeros(N,1)];
results.ebcmax = [ps zeros(N,1)];

box1 = (116:119)-62;
box2 = (148:155)-62;
ebc = (157:166)-62;

% Find minimum in box regions and maximum in EBC region
for i=1:N
    results.box1min(i,3) = min(bt(box1,i));
    results.box2min(i,3) = min(bt(box2,i));
    results.ebcmax(i,3)= max(bt(ebc,i));
end

% Associate values with their number and sort them, such that first entry
% is lowest/highest
results.box1Rank = sortrows([results.box1min(:,3) (1:N)'],1);
results.box2Rank = sortrows([results.box2min(:,3) (1:N)'],1);
results.ebcRank = sortrows([results.ebcmax(:,3) (1:N)'],-1);

% Put rankings back in order of pair swap indices
z1 = [results.box1Rank(:,2) (1:N)'];
z1 = sortrows(z1,1);
z2 = [results.box2Rank(:,2) (1:N)'];
z2 = sortrows(z2,1);
z3 = [results.ebcRank(:,2) (1:N)'];
z3 = sortrows(z3,1);

% Find score as sum of rankings
results.sumScore = sum([z1(:,2) z2(:,2) z3(:,2)]')';

% Now rank so that the top score shows first
results.sumScoreRanked = sortrows([results.sumScore (1:N)'],1);

for i=1:N
    results.sumScoreRanked(i,3:4) = ps(results.sumScoreRanked(i,2),:);
end

% 'results' will display an Nx4 matrix.  The first column is the score, the
% second column is the index of the pair swap, the third and fourth
% columnsn are the pair swap indices.
end


