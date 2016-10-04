function [results] = runSMCINModel(n1Rates, n4rates, xgrid, nTrials, nSteps)

% n1Rate = 0.2
% n4rates = 0.2 * ones(6*23,1);
N = 23;
%results.possStates = 1:5*N;
results.states = zeros(nTrials,nSteps);
results.states(:,1) = 2*N*ones(nTrials,1);

for i=1:nTrials
    for j=2:nSteps
        x = rand(1);
        
        % You can gain a chromosome, lose a chromosome, jump to 4N, or stay
        % the same
        if isnan(results.states(i,j-1))
            results.states(i,j) = NaN;
        elseif x < n1Rates(results.states(i,j-1))
            results.states(i,j) = results.states(i,j-1)+1;
        elseif x < 2*n1Rates(results.states(i,j-1))
            results.states(i,j) = results.states(i,j-1)-1;
        elseif x < (2*n1Rates(results.states(i,j-1)) + n4rates(results.states(i,j-1)))
            results.states(i,j) = 4*N;
        else
            results.states(i,j) = results.states(i,j-1);
        end
        if results.states(i,j) > (max(xgrid)-1)
            results.states(i,j) = NaN;
        elseif results.states(i,j) < 2
            results.states(i,j) = NaN;
        end
    end
end

end