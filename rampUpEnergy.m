function [results] = rampUpEnergy()

orgSeq = 'MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG';
rseq1 = orgSeq(randperm(length(orgSeq)));
rseq2 = orgSeq(randperm(length(orgSeq)));

increment = 1:10;
numTimes = 20;
results.SelfSelfScore = zeros(length(increment),1);
results.ABScore = zeros(length(increment),1);
[~,~,~,emin,~,~,~,~] = pholderAgg(rseq1,rseq1);
results.SelfSelfEmin = emin;
[~,~,~,emin,~,~,~,~] = pholderAgg(rseq1,rseq2);
results.ABEmin = emin;

for i=1:length(increment)
   disp(i);
   [score,~,~] = pholderVariance2Seq(rseq1,rseq1,numTimes,increment(i)); 
   results.SelfSelfScore(i) = score;
   
   [score,~,~] = pholderVariance2Seq(rseq1,rseq2,numTimes,increment(i)); 
   results.ABScore(i) = score;
end
end