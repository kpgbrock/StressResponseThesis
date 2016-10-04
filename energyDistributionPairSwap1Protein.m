function [results] = energyDistributionPairSwap1Protein()

N = 50;
nEnsemble = 200;

% Sperm whale myoglobin
orgSeq = 'MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG';

results.eminActualSelf = zeros(N,1);
results.eminActualAB = zeros(N,1);
results.eminSelf = zeros(nEnsemble,N);
results.eminAB = zeros(nEnsemble,N);
results.seq1 = cell(N,1);
results.seq2 = cell(N,1);
results.info = 'Sperm whale myoglobin permutations ensemble study.  Columns are ensembles, each row is a different permutation set.';


for i=1:N
    disp(i);
    seq1 = orgSeq(randperm(length(orgSeq)));
    seq2 = orgSeq(randperm(length(orgSeq)));
    results.seq1{i} = seq1;
    results.seq2{i} = seq2;
    
    % Get self-self energies
    [~,~,~,emin,~,~,~,~] = pholderAgg(seq1,seq1);
    results.eminActualSelf(i) = emin;       
    
    for k=1:nEnsemble
        seqSwap = swapMutantOneProtein(seq1);
        [~,~,~,emin,~,~,~,~] = pholderAgg(seqSwap,seqSwap);
        results.eminSelf(k,i) = emin;
    end
    
    [~,~,~,emin,~,~,~,~] = pholderAgg(seq1,seq2);
    results.eminActualAB(i) = emin; 
       
    for k=1:nEnsemble
        seqSwap1 = swapMutantOneProtein(seq1);
        seqSwap2 = swapMutantOneProtein(seq2);
        [~,~,~,emin,~,~,~,~] = pholderAgg(seqSwap1,seqSwap2);
        results.eminAB(k,i) = emin;
    end
    
end
end

% Input:  two sequences plus number of structures to generate
function [seq] = swapMutantOneProtein(seq)

n = length(seq);
z1 = randi(n);
z2 = z1;

while z2 == z1
    z2 = randi(n);
end

temp = seq(z1);
seq(z1) = seq(z2);
seq(z2) = temp;
end