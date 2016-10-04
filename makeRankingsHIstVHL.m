vhlGetRankingsAndScoresUntrunc;
r = rankOrderPairSwapVHLUntrunc;

[y1,x1] = hist(r.sumScoreRanked(:,1),100);
bar(x1,y1,'r');
hold on;
%bar(results(:,2),800*ones(20,1),'b');
plot(results(:,2),775*ones(20,1),'yv','MarkerSize',15,'MarkerFaceColor','y','MarkerEdgeColor','k','LineWidth',2.5);

axis fill;
lw = 6;
axis([-1000 70000 -1*lw 800]);
get(gca);
set(gca,'LineWidth',lw);
set(gca,'XTick',[0 3*10^4 6*10^4],'XTickLabel',['' '' ''],'YTick',[0 700], 'YTickLabel',['' '']);
%set(gca,'layer','top');
set(gca,'box','off');