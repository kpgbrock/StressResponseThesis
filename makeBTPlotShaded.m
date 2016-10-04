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

baseLine = 0;
greyColor = [0.5 0.5 0.5];
plot(x,maxBT,'Color',greyColor);
hold on;
h1 = fill(x([1 1:end end]), [baseLine maxBT' baseLine], greyColor, 'EdgeColor', 'none');
plot(x,minBT,'Color',greyColor);                             
h2 = fill(x([1 1:end end]),[baseLine minBT' baseLine],'w','EdgeColor','none');

plot(x,vhlUP.btTrunc(:,1)/max(vhlUP.btTrunc(:,1)),'k-.','LineWidth',4);
plot(x,vhlUP.btTrunc(:,20)/max(vhlUP.btTrunc(:,20)),'k','LineWidth',4);
