function [] = calcSameEPerm()

inputfile = 'yeastCytoplasmFiltered2.txt';
np = 10;
prec = 0.04;

data = fastaread(inputfile);

N = length(data);
ABnew = nan(N,1);
permABnew = nan(N,np);

whichA = 56;
Aseq = data(whichA,1).Sequence;


for i=1:N
    disp(i);
    Bseq = data(i,1).Sequence;
    [~,~,~,eminA,~,~,~,~] = pholderAgg(Aseq,Bseq);
    ABnew(i) = eminA;
    
    ct = 1;
    while ct <= np
        pA = Aseq(randperm(length(Aseq)));
        [~,~,~,emin,~,~,~,~] = pholderAgg(pA,Bseq);
        if (emin <= (eminA+prec) ) && (emin >= (eminA-prec))
            permABnew(i,ct) = emin;
            ct = ct + 1;
        end
    end
end

save ABnew ABnew;
save permABnew permABnew;
end