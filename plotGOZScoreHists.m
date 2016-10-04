function [results] = plotGOZScoreHists()

close all;
load goRP;
load goIsing100;
load godata;

N = length(godata.pIdx);

results.rpZNeutral = zeros(N,1);
results.rpZAcidic = zeros(N,1);
results.rpZDE = zeros(N,1);

results.isZNeutral = zeros(N,1);
results.isZAcidic = zeros(N,1);
results.isZDE = zeros(N,1);

for i=1:100
    results.rpZNeutral(i) = calcZScore(goRP.histNeutral(:,i),goRP.actNeutral(i));
    results.rpZAcidic(i) = calcZScore(goRP.histAcidic(:,i),goRP.actAcidic(i));
    results.rpZDE(i) = calcZScore(goRP.histDE(:,i),goRP.actDE(i));
    
    results.isZNeutral(i) = calcZScore(goIsing100.histNeutral(:,i),goIsing100.actNeutral(i));
    results.isZAcidic(i) = calcZScore(goIsing100.histAcidic(:,i),goIsing100.actAcidic(i));
    results.isZDE(i) = calcZScore(goIsing100.histDE(:,i),goIsing100.actDE(i));
end

fs = 20;
lw = 4;
x = -8:1:8;

subplot(1,2,1);
[h1,x1] = hist(results.rpZNeutral,x);
[h2,x2] = hist(results.rpZAcidic,x);
[h3,x3] = hist(results.rpZDE,x);
[h4,x4] = hist(results.isZNeutral,x);
[h5,x5] = hist(results.isZAcidic,x);
[h6,x6] = hist(results.isZDE,x);

h1 = h1/sum(h1);
h2 = h2/sum(h2);
h3 = h3/sum(h3);
h4 = h4/sum(h4);
h5 = h5/sum(h5);
h6 = h6/sum(h6);

plot(x1,h1,'b',x2,h2,'r',x3,h3,'k:','LineWidth',lw);

xlabel('Z Score');
ylabel('Counts');
axis([min(x) max(x) 0 1]);
title('Random Cytosolic Proteins','FontSize',fs);
set(gca,'FontSize',fs,'LineWidth',lw);
box off;
legend('Neutral','Acidic','\Delta E');
legend boxoff;

subplot(1,2,2);
plot(x4,h4,'b',x5,h5,'r',x6,h6,'k:','LineWidth',lw);

xlabel('Z Score');
ylabel('Counts');
axis([min(x) max(x) 0 1]);
title('Randomized Sequences (Same AA Content)','FontSize',fs);
set(gca,'FontSize',fs,'LineWidth',lw);
box off;
legend('Neutral','Acidic','\Delta E');
legend boxoff;

% subplot(3,2,1);
% h = histogram(results.rpZNeutral,x);
% h.Normalization = 'countdensity';
% xlabel('Neutral Z Score');
% ylabel('Counts');
% axis([-5 5 0 70]);
% title('Random Cytosolic Proteins','FontSize',fs);
% set(gca,'FontSize',fs,'LineWidth',lw);
% box off;
% 
% 
% subplot(3,2,3);
% h = histogram(results.rpZAcidic,x);
% h.Normalization = 'countdensity';
% xlabel('Acidic Z Score');
% ylabel('Counts');
% axis([-5 5 0 70]);
% set(gca,'FontSize',fs,'LineWidth',lw);
% box off;
% 
% subplot(3,2,5);
% h = histogram(results.rpZDE,x);
% h.Normalization = 'countdensity';
% xlabel('\Delta E Z Score');
% ylabel('Counts');
% axis([-5 5 0 70]);
% set(gca,'FontSize',fs,'LineWidth',lw);
% box off;
% 
% subplot(3,2,2);
% h = histogram(results.isZNeutral,x);
% h.Normalization = 'countdensity';
% xlabel('Neutral Z Score');
% ylabel('Counts');
% axis([-5 5 0 70]);
% title('Randomized Sequences (Same AA Content)','FontSize',fs);
% set(gca,'FontSize',fs,'LineWidth',lw);
% box off;
% 
% subplot(3,2,4);
% h = histogram(results.isZAcidic,x);
% h.Normalization = 'countdensity';
% xlabel('Acidic Z Score');
% ylabel('Counts');
% axis([-5 5 0 70]);
% set(gca,'FontSize',fs,'LineWidth',lw);
% box off;
% 
% 
% subplot(3,2,6);
% h = histogram(results.isZDE,x);
% h.Normalization = 'countdensity';
% xlabel('\Delta E Z Score');
% ylabel('Counts');
% axis([-5 5 0 70]);
% set(gca,'FontSize',fs,'LineWidth',lw);
% box off;

end

function [z] = calcZScore(x,realx)

z = (realx-mean(x))/std(x);

end