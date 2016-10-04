function [] = ptComparisonOriginal()

reps = 100;
FASTA = 'yeastCytoplasmFiltered2.txt';
rawdata = fastaread(FASTA);

data = zeros(reps*length(rawdata),1);

for j=1:length(rawdata)
    seq = upper(rawdata(j).Sequence);
    dataOrig = calcPt(seq);
    
    numMut = floor(length(seq)*0.1);
    mutLocs = randperm(length(seq));
    mutLocs = mutLocs(1:numMut);
    
    aa = ['I','V','L','F','C','M','A','G','T','W','S','Y','P','H','E','Q','D','N','K','R'];
    
    i = 1;
    while i <= reps
        seqMut = seq;
        for k=1:numMut
            inZ = mutLocs(k);
            repZ = floor(rand(1)*length(aa))+1;
            while seqMut(inZ) == aa(repZ) 
                repZ = floor(rand(1)*length(aa))+1;
            end
            seqMut(inZ) = aa(repZ);
        end
        data((j-1)*reps+i) = calcPt(seqMut)-dataOrig;
        i = i+1;
        
    end
end

save data data;

end

function [d] = calcPt(seq)

[copt,V,~,~,~,~,~,KDseq] = pholder(seq);
d = sum(abs(V*copt)'.*KDseq)/length(seq);

end