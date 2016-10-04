function [pStrict,pAE] = pValueSamplingPointMutVHL(test,nTrials,n)

countAbove = 0;
countEqual = 0;

for i=1:nTrials
    
   x = randi(2840,n,1);
   while length(unique(x)) ~= n
       x = randi(2840,n,1);
   end
   
   x = sum(x);
   
   if x > test
       countAbove = countAbove + 1;
   elseif x == test
       countEqual = countEqual + 1;
   end
end
pStrict = countAbove/nTrials;
pAE = (countAbove + countEqual)/nTrials;
end