close all;
load upVHLPairSwapResults;

m = zeros(142,142);
pairs = upVHLPairSwapResults.pairSwap;
cCrys = upVHLPairSwapResults.corrCrystal;

for i=1:length(pairs)
    m(pairs(i,1),pairs(i,2)) = cCrys(i);
    m(pairs(i,2),pairs(i,1)) = cCrys(i);
end

pcolor(63:204,63:204,m);
shading interp;
colorbar;