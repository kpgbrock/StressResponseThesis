x = linspace(-0.003,0.003,10);
close all;
load yeastEssential;
load y2hdata;

rBoth = isingModelIntrxns(y2hdata.idxNonSelf(yeastEssential.bothEssential,:),7,5);
rEither = isingModelIntrxns(y2hdata.idxNonSelf(yeastEssential.xorEssential,:),7,5);
rNeither = isingModelIntrxns(y2hdata.idxNonSelf(yeastEssential.neitherEssential,:),7,5);
% disp(rBoth.eTotalNeutral/length(rBoth.eComponentsNeutral));
% disp(rBoth.eTotalAcidic/length(rBoth.eComponentsNeutral))
% disp(rBoth.dE/length(rBoth.eComponentsNeutral));
% disp(rEither.eTotalNeutral/length(rEither.eComponentsNeutral));
% disp(rEither.eTotalAcidic/length(rEither.eComponentsNeutral));
% disp(rEither.dE/length(rEither.eComponentsNeutral));
% disp(rNeither.eTotalNeutral/length(rNeither.eComponentsNeutral));
% disp(rNeither.eTotalAcidic/length(rNeither.eComponentsNeutral));
% disp(rNeither.dE/length(rNeither.eComponentsNeutral));
disp(std(rBoth.eComponentsNeutral));
disp(std(rBoth.eComponentsAcidic))
disp(std(rBoth.eComponentsAcidic-rBoth.eComponentsNeutral));
disp(std(rEither.eComponentsNeutral));
disp(std(rEither.eComponentsAcidic));
disp(std(rEither.eComponentsAcidic-rEither.eComponentsNeutral));
disp(std(rNeither.eComponentsNeutral));
disp(std(rNeither.eComponentsAcidic));
disp(std(rNeither.eComponentsAcidic-rNeither.eComponentsNeutral));

[h,p] = kstest2(rBoth.eComponentsNeutral,rNeither.eComponentsNeutral)
[h,p] = kstest2(rBoth.eComponentsAcidic,rNeither.eComponentsAcidic)
[h,p] = kstest2(rBoth.eComponentsAcidic - rBoth.eComponentsNeutral,rNeither.eComponentsAcidic - rNeither.eComponentsNeutral)


[y1,~] = hist(rBoth.eComponentsNeutral,x);
[y2,~] = hist(rBoth.eComponentsAcidic,x);
[y3,~] = hist(rBoth.eComponentsAcidic - rBoth.eComponentsNeutral,x);

[y4,~] = hist(rEither.eComponentsNeutral,x);
[y5,~] = hist(rEither.eComponentsAcidic,x);
[y6,~] = hist(rEither.eComponentsAcidic - rEither.eComponentsNeutral,x);

[y7,~] = hist(rNeither.eComponentsNeutral,x);
[y8,~] = hist(rNeither.eComponentsAcidic,x);
[y9,~] = hist(rNeither.eComponentsAcidic - rNeither.eComponentsNeutral,x);

subplot(3,1,1);
plot(x,y1,'b',x,y2,'ro-',x,y3,'k^:');
legend('Neutral','Acidic','dE');
title('Both Proteins Essential');

subplot(3,1,2);
plot(x,y4,'b',x,y5,'r',x,y6,'k');
legend('Neutral','Acidic','dE');
title('1 Protein Essential');

subplot(3,1,3);
plot(x,y7,'b',x,y8,'r',x,y9,'k');
legend('Neutral','Acidic','dE');
title('Neither Protein Essential');

figure(2);
subplot(3,1,1);
lw = 4;
fs = 20;
plot(x,y1/sum(y1),'b',x,y4/sum(y4),'ro--',x,y7/sum(y7),'k^:','LineWidth',lw);

xlabel('E_{Neutral}');
ylabel('Counts');
box off;
set(gca,'LineWidth',lw,'FontSize',fs);

subplot(3,1,2);
plot(x,y2/sum(y2),'b',x,y5/sum(y5),'ro--',x,y8/sum(y8),'k^:','LineWidth',lw);
%legend('Both','Either','Neither');
%title('Acidic');
xlabel('E_{Acidic}');
ylabel('Counts');
box off;
set(gca,'LineWidth',lw,'FontSize',fs);

subplot(3,1,3);
plot(x,y3/sum(y3),'b',x,y6/sum(y6),'ro--',x,y9/sum(y9),'k^:','LineWidth',lw);
%legend('Both','Either','Neither');
xlabel('\Delta E');
ylabel('Counts');
box off;
set(gca,'LineWidth',lw,'FontSize',fs);