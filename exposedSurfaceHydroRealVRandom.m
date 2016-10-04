function [results] = exposedSurfaceHydroRealVRandom()

addpath 'C:/Users/Kelly/Desktop/EnglandRotation';

load btresults;
load yCyto50;
btreal = btresults.btReal;
btrand = btresults.btRandom(:,3);

avgshre = zeros(length(btreal),1);
numaboveshre = zeros(length(btreal),1);
percentaboveshre = zeros(length(btreal),1);

avgshra = zeros(length(btreal),1);
numaboveshra = zeros(length(btreal),1);
percentaboveshra = zeros(length(btreal),1);

hydrothreshold = 0;
msr = 0.75;
for i=1:length(btreal)
    re = btreal{i};
    ra = btrand{i};
    
    re = re(:,1);
    ra = ra(:,1);
    
    re = re > msr;
    ra = ra > msr;
    
    re = calcKD(yCyto50.sequence{i}(re));
    ra = calcKD(yCyto50.sequence{i}(ra));
    
    
    avgshre(i) = mean(re);
    numaboveshre(i) = length(re(re > hydrothreshold));
    percentaboveshre(i) = length(re(re > hydrothreshold))/length(re);
    
    avgshra(i) = mean(ra);
    numaboveshra(i) = length(ra(ra > hydrothreshold));
    percentaboveshra(i) = length(ra(ra > hydrothreshold))/length(ra);
end

results.realAvg = avgshre;
results.realNumAbove = numaboveshre;
results.realPercentAbove = percentaboveshre;

results.randomAvg = avgshra;
results.randomNumAbove = numaboveshra;
results.randomPercentAbove = percentaboveshra;

subplot(2,3,1);
hist(results.realAvg);
title('Real - Avg Exposed Hydro');
subplot(2,3,2);
hist(results.realNumAbove);
title('Real - Number Above Hydro Threshold');
subplot(2,3,3);
hist(results.realPercentAbove);
title('Real - Percent Above Hydro Threshold');
subplot(2,3,4);
hist(results.randomAvg);
title('Random - Avg Exposed Hydro');
subplot(2,3,5);
hist(results.randomNumAbove);
title('Random - Number Above Hydro Threshold');
subplot(2,3,6);
hist(results.randomPercentAbove);
title('Random - Percent Above Hydro Threshold');
end