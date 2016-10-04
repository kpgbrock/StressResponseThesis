close all;
load vhlUP;
mutsNot19 = vhlUP.btTrunc(:,[2:19,21]);
for i=1:19
    mutsNot19(:,i) = mutsNot19(:,i)/max(mutsNot19(:,i));
end

minBT = zeros(142,1);
maxBT = zeros(142,1);

for i=1:142
    minBT(i) = min(mutsNot19(i,:));
    maxBT(i) = max(mutsNot19(i,:));
end

[dist21VCBC,pdbseq] = createCrystalBurialTraceChain('1VCB',vhlUP.orgSeq(63:204),'C');
dist21VCBC = dist21VCBC/max(dist21VCBC);

hold on;

chap1 = (116:119);
chap2 = (148:155);
ebc = (157:166);

yellow = [255 207 0]/255;
blue = [0 0 191]/255;
red = [255 16 0]/255;
green = [0 204 0]/255;
% h1 = area(chap1,ones(length(chap1),1));
% h2 = area(chap2,ones(length(chap2),1));
% h3 = area(ebc,ones(length(ebc),1));
fill([chap1 chap1(end:-1:1)],[zeros(1,length(chap1)), ones(1,length(chap1))],yellow,'EdgeColor','none');
fill([chap2 chap2(end:-1:1)],[zeros(1,length(chap2)), ones(1,length(chap2))],blue,'EdgeColor','none');
fill([ebc ebc(end:-1:1)],[zeros(1,length(ebc)), ones(1,length(ebc))],red,'EdgeColor','none');

x = 63:204;

baseLine = 0;
greyColor = [0.5 0.5 0.5];

fill([x x(end:-1:1)],[minBT' maxBT(end:-1:1)'],greyColor);
plot(x,maxBT,'Color',greyColor);
plot(x,minBT,'Color',greyColor);                             

lw = 15;
plot(x,dist21VCBC,'Color',green,'LineWidth',lw);
plot(x,vhlUP.btTrunc(:,1)/max(vhlUP.btTrunc(:,1)),'ko','LineWidth',lw);
plot(x,vhlUP.btTrunc(:,20)/max(vhlUP.btTrunc(:,20)),'k','LineWidth',lw);

axis fill;
get(gca);
set(gca,'LineWidth',lw);
set(gca,'XTick',[120 180],'YTick',[0.5 1],'XTickLabel',['' ''],'YTickLabel',['','']);





