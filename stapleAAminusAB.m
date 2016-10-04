function [results] = stapleAAminusAB()

% Get sequences
load yCyto100300;
yc = yCyto100300;

numStaples = 10;

N = length(yc.sequence);

% Get indices of sequences to use as 'seqB' and make sure that it is not
% equivalent to seqA
B = randperm(N);
z = B - (1:N);
while ~isempty(z(z==0))
    B = randperm(N);
    z = B - (1:N);
end

% Create space for results
results.indexB = B;
results.eminAA = zeros(N,1);
results.eminAB = zeros(N,1);
results.diff = zeros(N,1);

for i=1:N
    disp(i);
    seqA = yc.sequence{i};
    seqB = yc.sequence{B(i)};
    nA = length(seqA);
    nB = length(seqB);
    
    
    [~,~,~,eminAA,~,~,~,~] = pholderStapleAB(seqA,seqA,makePairs(nA,nA,numStaples));
    [~,~,~,eminAB,~,~,~,~] = pholderStapleAB(seqA,seqB,makePairs(nA,nB,numStaples));
    
    results.eminAA(i) = eminAA;
    results.eminAB(i) = eminAB;
    results.diff(i) = eminAA - eminAB;
end
end

function [pairs] = makePairs(nA,nB,numStaples)

pairs = zeros(numStaples,2);

for j=1:numStaples
    
    pairs(j,1) = randi(nA);
    pairs(j,2) = randi(nB);
    pairs(j,2) = pairs(j,2)+nA;
end
end