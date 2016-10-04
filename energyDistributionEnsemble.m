function [results] = energyDistributionEnsemble(N)

load yCyto100300;
yc = yCyto100300;
L = 30;

% Get indices for two for AB (L1 and L2) and L3 for self-self
L1 = randperm(length(yc.abundance));
L1 = L1(1:L);

L2 = randperm(length(yc.abundance));
L2 = L2(1:L);

L3 = randperm(length(yc.abundance));
L3 = L3(1:L);


results.idAB1 = L1;
results.idAB2 = L2;
results.idSelf = L3;

results.eminAB = zeros(N,L);
results.eminActualAB = zeros(L,1);

results.eminSelf = zeros(N,L);
results.eminActualSelf = zeros(L,1);

% results.ensembleAB = cell(L,1);
% results.ensembleSelf = cell(L,1);

z = 1;
for i=1:L
    
    disp(i);
    % Remove star at end
    seq1 = yc.sequence{L1(i)};
    seq1(end) = [];
    
    seq2 = yc.sequence{L2(i)};
    seq2(end) = [];
    
    seq3 = yc.sequence{L3(i)};
    seq3(end) = [];
    
    
    [~,~,~,emin,~,~,phi0,~] = pholderAgg(seq1,seq2);
    results.eminActualAB(i) = emin;
    
    
    
    % Generate and apply AB ensemble
    ensemble = genEnsemblePerm(seq1,seq2,N);
    
   % results.ensembleAB{i} = ensemble;
    seqHydro = calcKD([seq1 seq2]);    
    for k=1:N
        results.eminAB(k,i) = ensemble.stretchE(k) + phi0*seqHydro*ensemble.r2(:,k);
    end
    
    [~,~,~,emin,~,~,phi0,~] = pholderAgg(seq3,seq3);
    results.eminActualSelf(i) = emin;
    
    % Generate and apply self/self ensemble
    ensembleself = genEnsemblePerm(seq3,seq3,N);
 
    seqHydro = calcKD([seq3 seq3]);
    for k=1:N
        results.eminSelf(k,i) = ensembleself.stretchE(k) + phi0*seqHydro*ensembleself.r2(:,k);
    end
    
    
    
end
end

function [] = genEnsembleRand(seq1,seq2,N)

n1 = length(seq1);
n2 = length(seq2);


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