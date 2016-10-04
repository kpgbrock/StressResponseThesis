function [seq] = makePercentSwapMutation(seq,percentToMutate)

numToMutate = ceil(percentToMutate * length(seq));
locsToMutate = randperm(length(seq));
locsToMutate = locsToMutate(1:numToMutate);

mutateTo = locsToMutate(randperm(numToMutate));

for i=1:length(locsToMutate)
    s = locsToMutate(i);
    s2 = mutateTo(i);
    tempAA = seq(s);
    seq(s) = seq(s2);
    seq(s2) = tempAA;
    
end
end