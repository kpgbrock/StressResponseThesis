function [polymer] = calcPolymerization(numProteins,numConf,t)


polymer = nan(t,3);
proteinInd = randi(numProteins);
confInd = randi(numConf);

polymer(1,1) = proteinInd;
polymer(1,2) = confInd;
polymer(1,3) = numProteins*(proteinInd-1)+confInd;
seed = polymer(1,:);

M = makeMUniform(numProteins,numConf);
for i=1:t
    
    proteinInd = randi(numProteins);
    confInd = randi(numConf);    
    z = numProteins*(proteinInd-1)+confInd;
    
    p = M(seed(1,3),z);
    if rand(1) < p
        polymer(i,1) = proteinInd;
        polymer(i,2) = confInd;
        polymer(i,3) = z;
        
        seed = polymer(i,:);
    end
        
end



end

function [M] = makeMUniform(numProteins,numConf)

z = numProteins*numConf;
M = zeros(z,z);

for i=1:z
    for j=1:z
        M(i,j) = rand(1);
    end
end

end