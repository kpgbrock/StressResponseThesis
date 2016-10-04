function [results] = seqSimilarityProteome(psize)

load yCyto100300;
yc = yCyto100300;

N = length(yc.abundance);

seqs = cell(N,1);

% Break proteome into chunks
for i=1:N
    seq = yc.sequence{i};
    seq(end) = [];
    seqs{i} = breakIntoDomains(seq,psize);
end

results.real = cell(N*(N-1)/2,1);
results.realrand = cell(N*(N-1)/2,1);
results.rand = cell(N*(N-1)/2,1);

z = 1;
for i=1:N
    for j=(i+1):N
        
        disp(z);
        % Disregard last chunk to ensure all domains of equal size
        seq1 = seqs{i};
        seq1(end) = [];
        seq2 = seqs{j};
        seq2(end) = [];
        
        numPairs = length(seq1)*length(seq2);
        results.real{z} = zeros(numPairs,1);
        results.realrand{z} = zeros(numPairs,1);
        results.rand{z} = zeros(numPairs,1);
        
        for k=1:length(seq1)
            for l=1:length(seq2)
                
%                 results.real{z}((k-1)*length(seq2)+l) = corrHydroPlot(seq1{k},seq2{l});
%                 results.rand{z}((k-1)*length(seq2)+l) = corrHydroPlot(seq1{k}(randperm(length(seq1{k}))),seq2{l}(randperm(length(seq2{l}))));
%                 results.realrand{z}((k-1)*length(seq2)+l) = corrHydroPlot(seq1{k}(randperm(length(seq1{k}))),seq2{l});
                results.real{z}((k-1)*length(seq2)+l) = seqSimBLOSUM(seq1{k},seq2{l});
                results.rand{z}((k-1)*length(seq2)+l) = seqSimBLOSUM(seq1{k}(randperm(length(seq1{k}))),seq2{l}(randperm(length(seq2{l}))));
                results.realrand{z}((k-1)*length(seq2)+l) = seqSimBLOSUM(seq1{k}(randperm(length(seq1{k}))),seq2{l});
            end
        end
        z = z+1;
    end
end



end

function [corrH] = corrHydroPlot(seq1,seq2)

kd1 = calcKD(seq1);
kd2 = calcKD(seq2);
corrH = corr(kd1',kd2');

end

function [score] = seqSimBLOSUM(seq1,seq2)

% Transforms AA chars into ints that correspond to the order in the
% MATLAB-given BLOSUM matrix
seq1 = aa2int(seq1);
seq2 = aa2int(seq2);

bl = blosum(62);

score = 0;
for i=1:length(seq1)
    score = score + bl(seq1(i),seq2(i));
end

end