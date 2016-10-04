function [results] = phCalcPredicatedYeastHuman(qOrg)

load phycfull;
z5 = find(phycfull.pHRange == 5);
z7 = find(phycfull.pHRange == 7);

results.yOpp7 = find((phycfull.charge(z7,qOrg.intrxnsBasic(:,1)).*phycfull.charge(z7,qOrg.intrxnsBasic(:,2))) < 0);
results.qOpp7 = length(find((qOrg.charge(z7,qOrg.intrxnsBasHID(results.yOpp7,1)).*qOrg.charge(z7,qOrg.intrxnsBasHID(results.yOpp7,2))) < 0));


results.ySame7 = find((phycfull.charge(z7,qOrg.intrxnsBasic(:,1)).*phycfull.charge(z7,qOrg.intrxnsBasic(:,2))) > 0);
results.qSame7 = length(find((qOrg.charge(z7,qOrg.intrxnsBasHID(results.ySame7,1)).*qOrg.charge(z7,qOrg.intrxnsBasHID(results.ySame7,2))) > 0));


results.yOpp5 = find((phycfull.charge(z5,qOrg.intrxnsBasic(:,1)).*phycfull.charge(z5,qOrg.intrxnsBasic(:,2))) < 0);
results.qOpp5 = length(find((qOrg.charge(z5,qOrg.intrxnsBasHID(results.yOpp5,1)).*qOrg.charge(z5,qOrg.intrxnsBasHID(results.yOpp5,2))) < 0));

results.ySame5 = find((phycfull.charge(z5,qOrg.intrxnsBasic(:,1)).*phycfull.charge(z5,qOrg.intrxnsBasic(:,2))) > 0);
results.qSame5 = length(find((qOrg.charge(z5,qOrg.intrxnsBasHID(results.ySame5,1)).*qOrg.charge(z5,qOrg.intrxnsBasHID(results.ySame5,2))) > 0));
    
end




