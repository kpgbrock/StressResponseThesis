function [x] = smMinCIN()

load skydata;
x0 = 0.2*ones(length(skydata.xgrid),1);
x = fmincon(@(x) smMinCINFun(x,skydata.xgrid,skydata.pdfgrid),x0,[],[],[],[],zeros(length(x0),1),0.5*ones(length(x0),1));
end

function [y] = smMinCINFun(n1Rates,xgrid,pdfgrid)

r = runSMCINModel(n1Rates,normpdf(xgrid,46,2), xgrid,10000, 100);
a = hist(r.states(:,end),xgrid);
a = a/sum(a);
b = pdfgrid/sum(pdfgrid);
y = -1*corr(a',b');
end