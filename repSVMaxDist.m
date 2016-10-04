% Want to reproduce the structural variability calculations that produced
% Jeremy's initial rankings of VHL mutants.

function [results] = repSVMaxDist()

% Get original sequence
load vhlUP;
orgSeq = vhlUP.orgSeq;
n = length(vhlUP.orgSeq);

k = 1;

z = nchoosek(n,2);
results.eminGS = zeros(z,1);
results.divScores = zeros(z,1);
results.dotProd = zeros(z,1);
results.possMut = zeros(z,2);
results.btGS = zeros(n,z);
results.btNew = zeros(n,z);

options = optimset('Display','off');

% Go through every possible pair swap mutation
for i=2:n
    disp(i);
    for j=1:(i-1)
        %disp(k);
        mutSeq = orgSeq;
        mutSeq([i,j]) = mutSeq([j,i]);
        
        lb = zeros(n,1);
        [copt,V,~,emin,r2max,msr,~,~] = pholder(mutSeq);
        btGS = V*copt;
        [btNew,dotprod] = linprog(btGS,diag(ones(n,1)),ones(n,1)*r2max,ones(1,n),msr * n,lb,[],[],options);
        
        results.divScores(k) = sum(abs(btGS-btNew));
        results.eminGS(k) = emin;
        results.dotProd(k) = dotprod;
        results.possMut(k,:) = [i,j];
        results.btGS(:,k) = btGS;
        results.btNew(:,k) = btNew;
        
        k = k+1;
    end
end

end