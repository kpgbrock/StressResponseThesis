function [results] = addAAdECorrWindow(seq,wsize)

n = length(seq);
z = n - wsize + 1;

results.dE = nan(z,1);
results.corrWithPrev = nan(z,1);
results.E = nan(z,1);

[copt,V,~,emin,~,~,~,~] = pholder(seq(1:wsize));
btPrev = V*copt;
ePrev = emin;
results.E(1) = emin;

for i=2:(n-wsize+1)
    disp(i);
    [copt,V,~,emin,~,~,~,~] = pholder(seq(i:(i+wsize-1)));
    btCurr = V*copt;
    eCurr = emin;
    
    results.corrWithPrev(i) = corr(btPrev,btCurr);
    results.dE(i) = eCurr - ePrev;
    results.E(i) = eCurr;
    
    btPrev = btCurr;
    ePrev = eCurr;
end

end