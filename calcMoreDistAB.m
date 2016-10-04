function [] = calcMoreDistAB()

inputfile = 'yeastCytoplasmFiltered2.txt';
data = fastaread(inputfile);

N = length(data);

ABall = nan(N*N,2);
ABminmax = nan(N,4);
pABall = nan(N*N,2);
pABminmax = nan(N,4);
permAseq = cell(N,1);

AB = nan(N,1);
pAB = nan(N,1);

for j = 1:N
    disp(j);
    Aseq = data(j,1).Sequence;
    pA = Aseq(randperm(length(Aseq)));
    permAseq{j} = pA;

    for i=1:N
        
        Bseq = data(i,1).Sequence;
        [~,~,~,emin,~,~,~,~] = pholderAgg(Aseq,Bseq);
        AB(i) = emin;
        
        [~,~,~,emin,~,~,~,~] = pholderAgg(pA,Bseq);
        pAB(i) = emin;
    end
    
    % Store minimum and maximum energies 
    [ABmin,ABminloc] = min(AB);
    ABminmax(j,1) = ABminloc;
    ABminmax(j,2) = ABmin;
    [ABmax,ABmaxloc] = max(AB);
    ABminmax(j,3) = ABmaxloc;
    ABminmax(j,4) = ABmax;
    
    [pABmin,pABminloc] = min(pAB);
    pABminmax(j,1) = pABminloc;
    pABminmax(j,2) = pABmin;
    [pABmax,pABmaxloc] = max(pAB);
    pABminmax(j,3) = pABmaxloc;
    pABminmax(j,4) = pABmax;
    
    ABall(((j-1)*N+1):((j-1)*N+N)) = AB-ABmin;
    pABall(((j-1)*N+1):((j-1)*N+N)) = pAB-pABmin;
    
   
end

save ABall ABall;
save ABminmax ABminmax;
save pABall pABall;
save pABminmax pABminmax;
save permAseq permAseq;

end
