function [results] = proteinSelfCorrVAbundance(psize)

load yCyto100300;
yc = yCyto100300;

N = length(yc.abundance);

results.corrHydro = cell(N,1);
results.simBLOSUM = cell(N,1);
results.corr2Letter = cell(N,1);
results.abundance = yc.abundance;

for i=1:N
    % Get sequence and discard * at end
    disp(i);
    seq = yc.sequence{i};
    seq(end) = [];
    
    % Break each protein into chunks and discard end piece to make sure all
    % are psize in length
    domains = breakIntoDomains(seq,psize);
    domains(end) = [];
    
    z = length(domains);
    results.corrHydro{i} = zeros((z*(z-1)/2),1);
    results.simBLOSUM{i} = zeros((z*(z-1)/2),1);
    results.corr2Letter{i} = zeros((z*(z-1)/2),1);
    x = 1;
    
    for k=1:z
        for j=(k+1):z
            results.corrHydro{i}(x) = corrHydroPlot(domains{k},domains{j});
            results.simBLOSUM{i}(x) = seqSimBLOSUM(domains{k},domains{j});
            results.corr2Letter{i}(x) = corr2LetterAlphabet(domains{k},domains{j});
            x = x+1;
        end
    end
end

end

function [corrH] = corrHydroPlot(seq1,seq2)

kd1 = calcKD(seq1);
kd2 = calcKD(seq2);
corrH = corr(kd1',kd2');

end

function [corrH] = corr2LetterAlphabet(seq1,seq2)

kd1 = calcKD(seq1);
kd2 = calcKD(seq2);
kd1(kd1 > 0) = 1;
kd1(kd1 == 0) = 0;
kd1(kd1 < 0) = -1;

kd2(kd2 > 0) = 1;
kd2(kd2 == 0) = 0;
kd2(kd2 < 0) = -1;

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