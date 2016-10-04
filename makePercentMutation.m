function [seq] = makePercentMutation(seq,percentToMutate)

AA = 'ARNDCEQGHILKMFPSTWYV';
numToMutate = ceil(percentToMutate * length(seq));
locsToMutate = randperm(length(seq));
locsToMutate = locsToMutate(1:numToMutate);

for i=1:length(locsToMutate)
    s = locsToMutate(i);
    seq(s) = AA(randi(length(AA),1));
    
end
end