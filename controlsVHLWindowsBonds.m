function [results] = controlsVHLWindowsBonds(windowSizes)

load vhlUP;
results.windowSizes = windowSizes;

x = 63:204;
n = length(windowSizes);
lw = 3;

results.hydroWT = zeros(length(x),1);
results.hydro19 = results.hydroWT;
results.corrWTWindow = zeros(1,n);
results.corr19Window = zeros(1,n);

for i=1:length(x)
    results.hydroWT(i) = KD(vhlUP.orgSeq(x(i)));
    results.hydro19(i) = KD(vhlUP.mutSeq{19}(x(i)));
end

results.wtAvgHydro = zeros(length(x),n);
results.mutAvgHydro = zeros(length(x),n);

for i=1:n
    w = windowSizes(i);
    
    % choose median
    medW = ceil(w/2);
    
    for j=1:(length(x)-w+1)
%         results.wtAvgHydro(j:j+w-1,i) = mean(results.hydroWT(j:j+w-1))*ones(w,1);
%         results.mutAvgHydro(j:j+w-1,i) = mean(results.hydro19(j:j+w-1))*ones(w,1);
        results.wtAvgHydro(j+medW-1,i) = mean(results.hydroWT(j:j+w-1));
        results.mutAvgHydro(j+medW-1,i) = mean(results.hydro19(j:j+w-1));
        
        
    end
    
    if w ~= 1
        results.wtAvgHydro(1:medW-1,i) = mean(results.hydroWT(1:w))*ones(medW-1,1);
        results.wtAvgHydro(end-floor(w/2)+1:end,i) = mean(results.hydroWT(end-w+1:end))*ones(floor(w/2),1);
        
        results.mutAvgHydro(1:medW-1,i) = mean(results.hydro19(1:w))*ones(medW-1,1);
        results.mutAvgHydro(end-floor(w/2)+1:end,i) = mean(results.hydro19(end-w+1:end))*ones(floor(w/2),1);
    end
    
    results.corrWTWindow(i) = corr(results.wtAvgHydro(:,i),vhlUP.btTrunc(:,1)-mean(vhlUP.btTrunc(:,1)));
    results.corr19Window(i) = corr(results.mutAvgHydro(:,i),vhlUP.btTrunc(:,20)-mean(vhlUP.btTrunc(:,20)));
end

figure(1);
subplot(2,2,1);
plot(x,results.wtAvgHydro,'LineWidth',lw);
legend(num2str(windowSizes'));
xlabel('Residue Index');
ylabel('Average Hydrophobicity');
title('A.    Wildtype Sequence');
set(gca,'LineWidth',4,'FontSize',16);
box off;

subplot(2,2,3);
plot(x,vhlUP.btTrunc(:,1)-mean(vhlUP.btTrunc(:,1)),'k^-','LineWidth',lw);
xlabel('Residue Index');
ylabel('R^2/R^2_{max}');
title('C.     Wildtype Sequence');
set(gca,'LineWidth',4,'FontSize',16);
box off;

subplot(2,2,2);
plot(x,results.mutAvgHydro,'LineWidth',lw);
legend(num2str(windowSizes'));
xlabel('Residue Index');
ylabel('Average Hydrophobicity');
title('B.     VHL19');
set(gca,'LineWidth',4,'FontSize',16);
box off;

subplot(2,2,4);
plot(x,vhlUP.btTrunc(:,20)-mean(vhlUP.btTrunc(:,20)),'k^-','LineWidth',lw);
xlabel('Residue Index');
ylabel('R^2/R^2_{max}');
title('D.     VHL19');
set(gca,'LineWidth',4,'FontSize',16);
box off;

figure(2);
plot(windowSizes,results.corrWTWindow,'bo',windowSizes,results.corr19Window,'k*');
xlabel('Window Size');
ylabel('Correlation with Burial Trace');
set(gca,'LineWidth',6','FontSize',24);
box off;
legend('Wildtype Sequence','VHL19');
end

