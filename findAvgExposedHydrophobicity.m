function [avgExpHydro] = findAvgExposedHydrophobicity(bt,seq,cutoff)

kdseq = zeros(length(seq),1);
for i=1:length(seq)
    kdseq(i) = KD(seq(i));
end

bt = bt/max(bt);

avgExpHydro = mean(kdseq(bt > cutoff));
end