function [results] = findRankScoreAnyVariable(bt)


% Normalize the burial traces
N = length(bt);
for i=1:N
    bt(:,i) = bt(:,i)/max(bt(:,i));
end

lseq = length(bt(:,1));
if lseq ~= 142
    disp('Incorrect length!');
    exit;
end

box1 = (116:119)-62;
box2 = (148:155)-62;
ebc = (157:166)-62;


results.box1min = zeros(N,1);
results.box2min = zeros(N,1);
results.ebcmax = zeros(N,1);

% Find minimum in box regions and maximum in EBC region
for i=1:N
    results.box1min(i) = min(bt(box1,i));
    results.box2min(i) = min(bt(box2,i));
    results.ebcmax(i)= max(bt(ebc,i));
end

results.scores = results.box1min + results.box2min + (1-results.ebcmax);

end


