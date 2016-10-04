function [results] = redoInitialVHLCalc(eAbove)

load vhlUP;
z = nchoosek(1:length(vhlUP.orgSeq),2);
results.scores = zeros(length(z),3);
results.kTAbove = eAbove;
results.seq = cell(length(z),1);
results.btOpt = zeros(length(vhlUP.orgSeq),length(z));
results.btGS = zeros(length(vhlUP.orgSeq),length(z));

[copt,V,~,~,~,msr,~,~] = pholder(vhlUP.orgSeq);
btWT = V*copt;
btWT = btWT-mean(btWT);

k = 1.5;
M = 10000000; 
n = length(vhlUP.orgSeq); 
rho0 = 250/(4*pi*(4^3)/3); 
r2max = (3*n/rho0/4/pi)^(2/3);
phi0 = (250/n)^(2/3)/112.5; 
% bond = zeros(n,n);
% for i = 1:n,
%     for j = 1:n,
%         bond(i,j) = (-((j + 1) == i) - ((i + 1) == j) + 2* (i == j))*k+M;
%     end
% end
% bond(1,1) = bond(1,1) - k;
% bond(n,n) = bond(n,n) - k;
options = optimset('Display','off');
lb = zeros(n,1);


for i=1:length(z)
    
    if mod(i,100)==0
        disp(i)
    end
    mutSeq = vhlUP.orgSeq;
    mutSeq([z(i,1),z(i,2)]) = mutSeq([z(i,2),z(i,1)]);
    [copt2,V2,eps2,emin2,~,~,~,~] = pholder(mutSeq);
    results.seq{i} = mutSeq;
    
    
%     seq = zeros(1,n);
%     for j = 1:n,
%         seq(j) = KD(mutSeq(j));
%     end
%     H = bond + phi0 * diag(seq);
%     [psi, epsmat] = eig(H);
%     eps = diag(epsmat);
%     V = psi.^2;
    
    A = V2;
    b = ones(n,1)*r2max;
    Aeq = [ones(1,n); eps2'];
    beq =[msr * n; emin2+eAbove];
    
    btGSNew = V2*copt2;
    btGSNew = btGSNew - mean(btGSNew);
    %[x,fval,~,~,~] = linprog(btWT'*V2,A,b,Aeq,beq,lb,[],[],options);
    [x,fval,~,~,~] = linprog((btGSNew)'*V2,A,b,Aeq,beq,lb,[],[],options);
    results.scores(i,:) = [fval z(i,:)];
    results.btOpt(:,i) = V2*x;
    results.btGS(:,i) = V2*copt2;
    
end
end