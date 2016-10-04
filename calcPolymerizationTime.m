function [results] = calcPolymerizationTime(numConf,cutoffLength)


% polymer = nan(t,3);
% proteinInd = randi(numProteins);
% confInd = randi(numConf);
%
% polymer(1,1) = proteinInd;
% polymer(1,2) = confInd;
% polymer(1,3) = numProteins*(proteinInd-1)+confInd;
% seed = polymer(1,:);
%
% for i=1:t
%
%     proteinInd = randi(numProteins);
%     confInd = randi(numConf);
%     z = numProteins*(proteinInd-1)+confInd;
%
%     p = M(seed(1,3),z);
%     if rand(1) < p
%         polymer(i,1) = proteinInd;
%         polymer(i,2) = confInd;
%         polymer(i,3) = z;
%
%         seed = polymer(i,:);
%     end
%
% end

numProteins = length(numConf);
M = makeMUniform(numConf);
z = sum(numConf);

results.selfself = zeros(z,2);
results.AB = zeros(z,2*cutoffLength+1);
results.M = M;
results.numProteins = numProteins;
results.numConf = numConf;

%cutoffLength = 10;
for i=1:z
    results.selfself(i,1) = i;
    results.selfself(i,2) = getTimeToPoly(M(i,i)*ones(cutoffLength,1));
    
    ABprob = zeros(cutoffLength,1);
    for j=1:cutoffLength
        
        pA = randi(numProteins);
        pB = pA;
        while (pB == pA)
            pB = randi(numProteins);
        end
        
        pAconf = sum(numConf(1:(pA-1))) + randi(numConf(pA));
        pBconf = sum(numConf(1:(pB-1))) + randi(numConf(pA));
        
        results.AB(i,(j-1)*2+1) = pAconf;
        results.AB(i,(j-1)*2+2) = pBconf;
        ABprob(j) = M(pAconf,pBconf);
    end
    results.AB(i,end) = getTimeToPoly(ABprob);
end

end

function [M] = makeMUniform(numConf)

z = sum(numConf);
M = zeros(z,z);

for i=1:z
    for j=1:z
        M(i,j) = rand(1);
    end
end

end

function [t] = getTimeToPoly(probs)

t = 0;
z = length(probs);
for i=1:z
    doWeAdd = 1;
    
    while(doWeAdd > probs(i))
       doWeAdd = rand(1); 
       t = t+1;
    end
    
end
end