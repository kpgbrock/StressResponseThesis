function [results] = stapleWithHigherK()

orgSeq = 'MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG';
N = length(orgSeq);

numStaples = 10;
numTrials = 50;
numSeqs = 1;


results.eminAA = zeros(numSeqs,numTrials);
results.eminAB = zeros(numSeqs,numTrials);
results.eminA = zeros(numSeqs,1);
results.eminB = zeros(numSeqs,1);
results.AAdiff = zeros(numSeqs,1);
results.ABdiff = zeros(numSeqs,1);

for i=1:numSeqs
    disp(i);
    seqA = orgSeq(randperm(N));
    seqB = orgSeq(randperm(N));
    nA = length(seqA);
    nB = length(seqB);
        
    [~,~,~,eminA,~,~,~,~] = pholder(seqA);
    [~,~,~,eminB,~,~,~,~] = pholder(seqB);
    results.eminA(i) = eminA;
    results.eminB(i) = eminB;
    
    for j=1:numTrials
        [~,~,~,eminAA,~,~,~,~] = pholderStapleABnewK(seqA,seqA,makePairs(nA,nA,numStaples));
        [~,~,~,eminAB,~,~,~,~] = pholderStapleABnewK(seqA,seqB,makePairs(nA,nB,numStaples));
                
        results.eminAA(i,j) = eminAA;
        results.eminAB(i,j) = eminAB;
        
    end
    
    results.AAdiff(i) = min(results.eminAA(i,:))-2*eminA;
    results.ABdiff(i) = min(results.eminAB(i,:)) - eminA - eminB;
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