load isingResults10000Fixed;
x = linspace(-0.2,0.5,20);
lw = 4;
fs = 16;

subplot(3,1,1);
hist(isingResults10000Fixed.eTotalNeutral,x);
axis([-0.2 0.5 0 4000]);
hold on;
plot([isingResults10000Fixed.actETotalNeutral isingResults10000Fixed.actETotalNeutral], get(gca,'YLim'),'k');
box off;
xlabel('\Sigma E_{Neutral}');
ylabel('Counts');
set(gca,'LineWidth',lw,'FontSize',fs);
title('0.4% Less','FontSize',fs);

subplot(3,1,2);
hist(isingResults10000Fixed.eTotalAcidic,x);
axis([-0.2 0.5 0 4000]);
hold on;
plot([isingResults10000Fixed.actETotalAcidic isingResults10000Fixed.actETotalAcidic], get(gca,'YLim'),'k');
box off;
xlabel('\Sigma E_{Acidic}');
ylabel('Counts');
set(gca,'LineWidth',lw,'FontSize',fs);
title('0% Less','FontSize',fs);

subplot(3,1,3);
hist(isingResults10000Fixed.dEHist,x);
axis([-0.2 0.5 0 4000]);
hold on;
plot([isingResults10000Fixed.actdE isingResults10000Fixed.actdE], get(gca,'YLim'),'k');
box off;
xlabel('\Sigma \Delta E');
ylabel('Counts');
set(gca,'LineWidth',lw,'FontSize',fs);
title('11.7% Less','FontSize',fs);