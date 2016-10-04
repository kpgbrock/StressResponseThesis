function [seqs] = genNSequences(m)

seqs = cell(length(m),1);
for i=1:length(m)
    seqs{i} = '';
    for j=1:m(i)
        seqs{i}(j) = chooseRandAA();
    end
end
end

function [a] = chooseRandAA()
AA = 'ARNDCEQGHILKMFPSTWYV';
a = AA(randi(length(AA),1));
end