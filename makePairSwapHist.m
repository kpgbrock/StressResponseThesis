close all;
load vhlUP;
load upVHLPairSwapResults;

mutsNot19 = vhlUP.btTrunc(:,[2:19,21]);
for i=1:19
    mutsNot19(:,i) = mutsNot19(:,i)/max(mutsNot19(:,i));
end

nbins = 30;

indexMut19 = find((upVHLPairSwapResults.pairSwap(:,1) == (173-62)) & (upVHLPairSwapResults.pairSwap(:,2) == (201-62)));
hist(upVHLPairSwapResults.corrCrystal,nbins);

length(upVHLPairSwapResults.corrCrystal(upVHLPairSwapResults.corrCrystal > upVHLPairSwapResults.corrCrystal(indexMut19)))/length(upVHLPairSwapResults.corrCrystal)


lw = 15;
axis fill;
axis([-0.2 0.6 0 2200]);
get(gca);
set(gca,'LineWidth',lw);
set(gca,'XTick',[0 upVHLPairSwapResults.corrCrystal(indexMut19), 0.5],'XTickLabel',['' '' ''],'YTick',[0 2000], 'YTickLabel',['' '']);
set(gca,'layer','top');
set(gca,'box','off');





