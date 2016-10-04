function [results] = marketInd(Q,Fa,Fb,bA0,bB0,dA0,dB0,Mt0,calcMt,alpha,delD,delB,threshold,disc,stepMax)

% Want to count number of steps until reach steady state
steps = 0;

Mt = Mt0;

bA = bA0;
bB = bB0;
dA = dA0;
dB = dB0;

vecBA = zeros(stepMax,1);
vecBB = zeros(stepMax,1);
vecDA = zeros(stepMax,1);
vecDB = zeros(stepMax,1);

while(steps < stepMax)
    
    vecBA(steps+1) = bA;
    vecBB(steps+1) = bB;
    vecDA(steps+1) = dA;
    vecDB(steps+1) = dB;
    
   betaA = alpha/((1-0.95*(1-disc^(steps+1)))*Fa);
   betaB = alpha/((1-0.95*(1-disc^(steps+1)))*Fb);
   
   D_b_A = Mt*(Q*dA^e)/(Q*dA^e + Q*dB^e) - delB*bA;
   D_b_B = Mt*(Q*dB^e)/(Q*dA^e + Q*dB^e) - delB*bB;
   
   D_d_A = betaA*bA - delD*dA;
   D_d_B = betaB*bB - delD*dB;
   
   bA = bA+D_b_A;
   bB = bB+D_b_B;
   dA = dA + D_d_A;
   dB = dB + D_d_B;
   
   Mt = calcMt(Mt);
   steps = steps+1;
   
   if ((D_b_A < 10^threshold) && (D_b_B < 10^threshold) && (D_d_A < 10^threshold) && (D_d_B < 10^threshold) )
       vecBA = vecBA(1:steps);
       vecBB = vecBB(1:steps);
       vecDA = vecDA(1:steps);
       vecDB = vecDB(1:steps);
       break;
   end
   
end

results.bA = vecBA;
results.bB = vecBB;
results.dA = vecDA;
results.dB = vecDB;
results.numSteps = length(results.bA);

end