function [results] = addAAdECorr(seq,startPos,endPos)

n = length(seq);
z = endPos - startPos;

results.dE = zeros(z,1);
results.corrWithPrev = zeros(z,1);

[copt,V,~,emin,~,~,~,~] = pholder(seq(1:startPos));
btPrev = V*copt;
ePrev = emin;

for i=(startPos+1):endPos
    disp(i);
    [copt,V,~,emin,~,~,~,~] = pholder(seq(1:i));
    btCurr = V*copt;
    eCurr = emin;
    
    results.corrWithPrev(i-startPos) = corr(btPrev,btCurr(1:end-1));
    results.dE(i-startPos) = eCurr - ePrev;
    
    btPrev = btCurr;
    ePrev = eCurr;
end

end