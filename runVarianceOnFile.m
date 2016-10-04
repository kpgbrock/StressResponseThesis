% Kelly Brock
% Do pholderVariance for a list of items

function [final] = runVarianceOnFile(infile)

%infile = 'C:\Users\Kelly\Desktop\MMPS\PurkinjeCells.fasta';
data1 = fastaread(infile);

scoredata = nan(length(data1),1);
resultsdata = cell(length(data1),1);
efdata = cell(length(data1),1);

for i=1:length(data1)
   
   seq = data1(i,1).Sequence;
   if ((length(seq) >= 100) && (length(seq) <= 300))
       disp(i);
       [score,results,ef] = pholderVariance(seq,100);
       disp(score);
       scoredata(i) = score;
       resultsdata{i} = results;
       efdata{i} = ef;
   end
end

final.scores = scoredata;
final.bt = resultsdata;
final.ef = efdata;
end