function [pAbove,pAE] = pValuePointMutationVHL(test)
N = 2840;

countAbove = 0;
countEqual = 0;
countTotal = 0;
for i=1:N
    
    for j=[1:(i-1) (i+1):N]
        klist = 1:N;
        klist([i,j]) = [];
        j
        for k=klist
            
            llist = 1:N;
            llist([i,j,k]) = [];
            
            for l=llist
                
                s = i + j + k + l;
                if s == test
                    countEqual = countEqual+1;
                elseif s > test
                    countAbove = countAbove + 1;
                end
                countTotal = countTotal + 1;
            end
        end
    end
end

pAbove = countAbove/countTotal;
pAE = (countAbove + countEqual)/countTotal;
end
                
                