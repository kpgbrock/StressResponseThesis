function [results] = scatterplotABHydroEContributions(N)

orgSeq = 'MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG';
L = length(orgSeq);

results = zeros(N,2);
for i=1:N
   seq1 = orgSeq(randperm(length(orgSeq)));
   seq2 = orgSeq(randperm(length(orgSeq)));
   
   % Self/nonself
    [copt,V,~,~,~,~,phi0,seq] = pholderAgg(seq1,seq2);
    bt = V*copt;
    results(i,1) = phi0*seq(1:L)'*bt(1:L);
    results(i,2) = phi0*seq((L+1):end)'*bt((L+1):end);
    
end
end