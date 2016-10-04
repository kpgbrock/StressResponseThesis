function [results] = checkAllosteryByChangingHydro(seq,deltaHall)

[copt,V,~,eminOrg,~,~,~,~] = pholder(seq);
btOrg = V*copt;
btOrg = btOrg/max(btOrg);

N = length(seq);
%deltaHall = [-10,-5,-1,1,5,10];

results.delBT = zeros(N,N,length(deltaHall));
results.delE = zeros(N,length(deltaHall));
results.deltaHall = deltaHall;
results.seq = seq;

for i=1:N
    
    kdOrg = KD(seq(i));
    disp(i);
    for j=1:length(deltaHall)
        
        deltaH = kdOrg+deltaHall(j);
        %deltaH = deltaHall(j);
        [copt,V,~,emin,~,~,~,~] = pholderPinAAToSurface(seq,i,deltaH);
        btNew = V*copt;
        btNew = btNew/max(btNew);
        
        results.delBT(i,:,j) = btNew-btOrg;
        results.delE(i,j) = emin-eminOrg;
    end
end
end