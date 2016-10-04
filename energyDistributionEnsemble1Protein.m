function [results] = energyDistributionEnsemble1Protein()

N = 100;
nEnsemble = 1000;

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
    
    ensembleSelf = genEnsemblePerm(seq1,seq1,nEnsemble);
    ensembleRP1RP2 = genEnsemblePerm(seq1,seq2,nEnsemble);
    
    % Get self-self energies
    [~,~,~,emin,~,~,phi0,~] = pholderAgg(seq1,seq1);
    results.eminActualSelf(i) = emin;     
    seqHydro = calcKD([seq1 seq1])';  
    
    for k=1:nEnsemble
        results.eminSelf(k,i) = ensembleSelf.stretchE(k) + phi0*seqHydro*ensembleSelf.r2(:,k);
    end
    
    [~,~,~,emin,~,~,phi0,~] = pholderAgg(seq1,seq2);
    results.eminActualAB(i) = emin; 
    
    seqHydro = calcKD([seq1 seq2])';    
    for k=1:nEnsemble
        results.eminAB(k,i) = ensembleRP1RP2.stretchE(k) + phi0*seqHydro*ensembleRP1RP2.r2(:,k);
    end
    
end
end

% Input:  two sequences plus number of structures to generate
function [results] = genEnsemblePerm(seq1,seq2,N)

n1 = length(seq1);
n2 = length(seq2);

results.r2 = zeros(length(seq1)+length(seq2),N);
results.stretchE = zeros(N,1);

for i=1:N
    rseq1 = seq1(randperm(n1));
    rseq2 = seq2(randperm(n2));
    
    [copt,V,~,emin,~,~,phi0,seq] = pholderAgg(rseq1,rseq2);
    results.r2(:,i) = V*copt;
    results.stretchE(i) = emin - phi0*seq'*results.r2(:,i);
end
end