function [fcharge,pI] = findPICurveWikiNoEnds(seq,pHRange)



%pHRange = 1:0.01:14;



seq = upper(seq);

[~,resPos,resNeg] = findpKaWiki('');



N = length(seq);



fcharge = zeros(length(pHRange),1);



nPos = length(resPos);

resInt = strcat(resPos,resNeg);



% Go through the set of possible ionized side chains

for i=1:length(resInt)

   charInt = resInt(i);

   resCount = length(seq(seq==charInt));

   

   

   % Decide whether amino acid can be protonated or deprotonated

   neg = false;

   if i > nPos

       neg = true;

   end

   

   % Calculate charge accordingly

   if resCount > 0       

       fcharge = fcharge + calcFValues(resCount,pHRange,charInt,neg);

   end      

end





% Add in pK values for N and C terminus
% 
% fcharge = calcFValues(1,pHRange,strcat(seq(1),'N'),false) + fcharge;
% 
% fcharge = calcFValues(1,pHRange,strcat(seq(end),'C'),true) + fcharge;



% Find where sign changes for pI calculation (if needed)

indZero = find(fcharge(1:end-1).*fcharge(2:end) < 0);

pI = (pHRange(indZero) + pHRange(indZero+1))/2;



end



function [fcharge] = calcFValues(resCount,pHRange,charInt,neg)

   fcharge = (1./(1+10.^(pHRange - findpKaWiki(charInt))))';

   if neg

       fcharge = fcharge -1;

   end

   

   fcharge = fcharge*resCount;

end