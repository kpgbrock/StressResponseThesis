function [] = calcDistAB()

inputfile = 'yeastCytoplasmFiltered2.txt';
data = fastaread(inputfile);

N = length(data);

ABmin = nan(N,1);
ABmed = nan(N,1);
pABmin = nan(N,1);
pABmed = nan(N,1);

AB = nan(N,1);
pAB = nan(N,1);

for j = 1:N
    disp(j);
    Aseq = data(j,1).Sequence;
    pA = Aseq(randperm(length(Aseq)));

    for i=1:N
        
        Bseq = data(i,1).Sequence;
        [~,~,~,emin,~,~,~,~] = pholderAgg(Aseq,Bseq);
        AB(i) = emin;
        
        [~,~,~,emin,~,~,~,~] = pholderAgg(pA,Bseq);
        pAB(i) = emin;
    end
    
    ABmin(j) = min(AB);
    ABmed(j) = median(AB);
    pABmin(j) = min(pAB);
    pABmed(j) = median(pAB);
    
   
end

save ABmin ABmin;
save ABmed ABmed;
save pABmin pABmin;
ssave pABmed pABmed;
end
