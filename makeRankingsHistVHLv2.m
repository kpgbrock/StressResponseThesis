close all;
load upVHLPairSwapAll;
r = findRankScoreAnyVariable(upVHLPairSwapAll.bt);

load vhlUP;
z = vhlUP.pairSwaps;
locsOfMutants = zeros(20,1);
for i=1:20
    if z(i,1) > z(i,2)
        z(i,:) = [z(i,2) z(i,1)];
    end
    locsOfMutants(i) = find((upVHLPairSwapAll.pairSwap(:,1) == z(i,1)) & (upVHLPairSwapAll.pairSwap(:,2) == z(i,2)));
end


[y1,x1] = hist(r.scores,100);
hold on;
h = bar(x1,y1,'r');
ylimits = get(gca,'YLim');

for i=1:20
   
    plot([r.scores(locsOfMutants(i)) r.scores(locsOfMutants(i))],ylimits,'--','color',[0.5 0.5 0.5],'LineWidth',4);
 
end

uistack(h,'top');

%bar(results(:,2),800*ones(20,1),'b');
%plot(results(:,2),775*ones(20,1),'yv','MarkerSize',15,'MarkerFaceColor','y','MarkerEdgeColor','k','LineWidth',2.5);

axis fill;
lw = 6;
axis([0.4 2 0 7000]);
get(gca);
set(gca,'LineWidth',lw);
set(gca,'XTick',[0.5 1 1.5 2],'XTickLabel',['' '' '' ''],'YTick',[0 7000], 'YTickLabel',['' '']);
%set(gca,'XTick',[0.5 1 1.5 2],'YTick',[0 7000]);
%set(gca,'layer','top');
set(gca,'TickLength',[0.03 0.03]);
set(gca,'box','off');