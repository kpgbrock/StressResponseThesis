function [score,results,ef] = pholderVariance2Seq(seq1,seq2,numTimes,increment)

[copt,V,eps,emin,r2max,msr,~,~] = pholderAgg(seq1,seq2);
results = zeros(length(seq1)+length(seq2),numTimes);

results(:,1) = V*copt;
results(:,1) = results(:,1)/max(results(:,1));

%increment = 2;
k = 1.5; 
M = 10000000; 
n = length(seq1)+length(seq2); 

lb = zeros(n,1); 
A = V;
b = ones(n,1)*r2max;
Aeq = [ones(1,n); eps'];
beq =[msr * n; emin+increment];

ef = [];
for i=2:numTimes
    fNew = randn(n,1);
    fNew = (fNew./max(fNew))+eps;
    options = optimset('Display','off');
    [copt2,emin2,exitflag] = linprog(fNew,A,b,Aeq,beq,lb,[],[],options);
    if exitflag == 1
        results(:,i) = V*copt2;
        results(:,i) = results(:,i)/max(results(:,i));
    else
        ef = [ef; i];
    end    
end

results(:,ef) = [];

score = sum(std(results').^2)./n;

end