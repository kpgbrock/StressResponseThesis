% Want to reproduce the structural variability calculations that produced
% Jeremy's initial rankings of VHL mutants.

function [results] = repSVCalcsJE(nSV,increment)

% Get original sequence
load vhlUP;
orgSeq = vhlUP.orgSeq;

k = 1;

z = nchoosek(213,2);
results.divScores = zeros(z,1);
results.divScores = zeros(z,1);
results.possMuts = zeros(z,2);

% Go through every possible pair swap mutation
for i=2:213
    for j=1:(i-1)
        %disp(k);
        mutSeq = orgSeq;
        mutSeq([i,j]) = mutSeq([j,i]);
        
        [score,bts,~] = pholderVarianceInc(mutSeq,nSV,increment);
        
        divSc = zeros(nchoosek(nSV,2),1);
        n = 1;
        for l = 2:nSV
            for m = 1:(l-1)
                divSc(n) = sum(abs(bts(:,l)-bts(:,m)));
                n = n+1;
            end
        end
        
        results.divScores(k) = max(divSc);
        results.strScores(k) = score;
        k = k+1;
    end
end

end