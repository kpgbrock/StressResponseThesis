function [results] = stapleAAminusAB1Protein()

orgSeq = 'MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG';
N = length(orgSeq);

numStaples = 10;
numTrials = 500;

results.eminAA = zeros(N,1);
results.eminAB = zeros(N,1);
results.diff = zeros(N,1);

for i=1:numTrials
    disp(i);
    seqA = orgSeq(randperm(N));
    seqB = orgSeq(randperm(N));
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