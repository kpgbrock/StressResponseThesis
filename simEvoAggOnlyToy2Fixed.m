function [results] = simEvoAggOnlyToy2Fixed(seqlength,gens,percentToMutate)

load yCyto50;
yc = yCyto50;
N = length(yc.id);
for i=1:N
    tempseq = yc.sequence{i}(1:end-1);
    tempseq = tempseq(randperm(length(tempseq)));
    yc.sequence{i} = tempseq(1:seqlength);
end
seqs = yc.sequence;
[ap,selfAP] = calcAggProteomeNoSAH(seqs);
threshold = min([ap;selfAP]);

apFull = zeros(gens+1,N*(N-1)/2);
sapFull = zeros(gens+1,N);

avgAP = zeros(gens+1,1);
stdAP = zeros(gens+1,1);
maxAP = zeros(gens+1,1);
avgSAP = zeros(gens+1,1);
stdSAP = zeros(gens+1,1);
maxSAP = zeros(gens+1,1);

avgAP(1) = mean(ap);
stdAP(1) = std(ap);
maxAP(1) = max(ap);
avgSAP(1) = mean(selfAP);
stdSAP(1) = std(selfAP);
maxSAP(1) = max(selfAP);
apFull(1,:) = ap;
sapFull(1,:) = selfAP;

seqsOld = seqs;
apOld = ap;
selfAPOld = selfAP;

avgAPOld = avgAP(1);
stdAPOld = stdAP(1);
maxAPOld = maxAP(1);
avgSAPOld = avgSAP(1);
stdSAPOld = stdSAP(1);
maxSAPOld = maxSAP(1);

for i=1:gens
    disp(i);
    seqGen = swapMutateAllProteinsPercent(seqsOld,percentToMutate);
    [apGen,selfAPGen] = calcAggProteomeNoSAH(seqGen);
    newmax = min([apGen;selfAPGen]);
    if newmax > threshold        
        apOld = apGen;
        selfAPOld = selfAPGen;
        
        
        avgAPOld = mean(apGen);
        stdAPOld = std(apGen);
        maxAPOld = max(apGen);
        avgSAPOld = mean(selfAPGen);
        stdSAPOld = std(selfAPGen);
        maxSAPOld = max(selfAPGen);
        
        
    
        
    end
    avgAP(i+1) = avgAPOld;
    stdAP(i+1) = stdAPOld;
    maxAP(i+1) = maxAPOld;
    avgSAP(i+1) = avgSAPOld;
    stdSAP(i+1) = stdSAPOld;
    maxSAP(i+1) = maxSAPOld;
    
    apFull(i+1,:) = apOld';
    sapFull(i+1,:) = selfAPOld';
end

results.avgAP = avgAP;
results.stdAP = stdAP;
results.maxAP = maxAP;
results.avgSAP = avgSAP;
results.stdSAP = stdSAP;
results.maxSAP = maxSAP;

results.apFull = apFull;
results.sapFull = sapFull;
results.endSeqs = seqsOld;
results.firstSeqs = seqs;

end