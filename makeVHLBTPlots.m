close all;
figure(1);
load vhl;
x = 63:204;
box1x = 116:119;
box1=box1x-(x(1)-1);
box2x = 148:155;
box2=box2x-(x(1)-1);
ebcx = 157:166;
ebc = ebcx-(x(1)-1);
[dist21VCBC,pdbseq] = createCrystalBurialTraceChain('1VCB',vhl.orgSeq(x),'C');
hold on;

for i=2:19
    bt = vhl.btTrunc(:,i)/max(vhl.btTrunc(:,i));
    plot(x,bt,'k:');
end

i=21;
bt = vhl.btTrunc(:,i)/max(vhl.btTrunc(:,i));
plot(x,bt,'k:');

lw = 2;
lw2 = 3;
bt = vhl.btTrunc(:,1)/max(vhl.btTrunc(:,1));
plot(x,bt,'c','LineWidth',lw);
% plot(box1x,bt(box1),'cs-','LineWidth',lw2);
% plot(box2x,bt(box2),'co-','LineWidth',lw2);
% plot(ebcx,bt(ebc),'c^-','LineWidth',lw2);

bt = dist21VCBC/max(dist21VCBC);
plot(x,bt,'b','LineWidth',lw);
% plot(box1x,bt(box1),'bs-','LineWidth',lw2);
% plot(box2x,bt(box2),'bo-','LineWidth',lw2);
% plot(ebcx,bt(ebc),'b^-','LineWidth',lw2);

bt = vhl.btTrunc(:,20)/max(vhl.btTrunc(:,20));
plot(x,bt,'r','LineWidth',lw);
% plot(box1x,bt(box1),'rs-','LineWidth',lw2);
% plot(box2x,bt(box2),'ro-','LineWidth',lw2);
% plot(ebcx,bt(ebc),'r^-','LineWidth',lw2);

