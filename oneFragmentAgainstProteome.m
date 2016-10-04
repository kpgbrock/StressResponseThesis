function [results] = oneFragmentAgainstProteome(psize,seqIndex,domainBase)

load yCyto100300;
yc = yCyto100300;
n = length(yc.abundance);

% If an input fragment isn't given
if narg == 1
    
    % Choose one random fragment
    
    seqIndex = randi(n);
    seqBase = yc.sequence{seqIndex};
    seqBase(end) = [];
    seqBase = breakIntoDomains(seqBase,psize);
    seqBase(end) = [];
    domainIndex = randi(length(seqBase));
    domainBase = seqBase{domainIndex};
end

results.abundanceTemplate = yc.abundance(seqIndex);


yc.sequence(seqIndex) = [];
yc.abundance(seqIndex) = [];

n = n-1;

% Uncomment to randomly shuffle abundances
% yc.abundance = yc.abundance(randperm(length(yc.abundance)));

% Break remaining proteome into chunks
tally = 0;
seqs = cell(n,1);
for i=1:n
    
    seq = yc.sequence{i};
    seq(end) = [];
    
    % Uncomment to do for random permutations
    % seq = seq(randperm(length(seq)));
    
    seqs{i} = breakIntoDomains(seq,psize);
    seqs{i}(end) = [];
    tally = tally+yc.abundance(i)*length(seqs{i});
    
end

results.abundanceComparison = yc.abundance;
results.simCorrHydroDistro = zeros(tally,1);
results.simBlosumDistro = zeros(tally,1);
results.seqTemplate = seqIndex;
results.psize = psize;
results.domainTemplate = domainBase;
results.tally = tally;

z = 1;
for i=1:n  
    
    for j=1:length(seqs{i})
        blosum = seqSimBLOSUM(domainBase,seqs{i}{j});
        hydro = corrHydroPlot(domainBase,seqs{i}{j});
        for k=1:yc.abundance(i)
            results.simCorrHydroDistro(z) = hydro;
            results.simBlosumDistro(z) = blosum;
            z = z+1;
        end               
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