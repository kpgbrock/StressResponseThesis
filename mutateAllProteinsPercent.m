function [seqsMut] = mutateAllProteinsPercent(seqs,percentToMutate)

N = length(seqs);
seqsMut = cell(N,1);

for i=1:N
    seqsMut{i} = makePercentMutation(seqs{i},percentToMutate);
end
end