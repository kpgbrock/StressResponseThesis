function [ap,selfAP] = calcAggProteomeNoSAH(seqs)

N = length(seqs);
numAP = N*(N-1)/2;

ap = zeros(numAP,1);
selfAP = zeros(N,1);

z = 1;
for i=1:N
    for j=1:(i-1)
        [~,~,~,emin,~,~,~,~] = pholderAgg(seqs{i},seqs{j});
        ap(z) = emin;
        z = z+1;
    end
end

for i=1:N
    [~,~,~,emin,~,~,~,~] = pholderAgg(seqs{i},seqs{i});
    selfAP(i) = emin;
end
end