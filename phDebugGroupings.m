clear all;
load phIntns;

zOrg = phIntns.modChargeGroups;
zNew = phIntns.modChargeGroups;
n = zNew(:,5);
n = n(randperm(length(n)));
zNew(:,5) = n;

z1Org = zOrg(zOrg(:,5)==1,:);
z2Org = zOrg(zOrg(:,5)==2,:);
z3Org = zOrg(zOrg(:,5)==3,:);
z4Org = zOrg(zOrg(:,5)==4,:);

z1New = zNew(zNew(:,5)==1,:);
z2New = zNew(zNew(:,5)==2,:);
z3New = zNew(zNew(:,5)==3,:);
z4New = zNew(zNew(:,5)==4,:);

uOrg1 = [z1Org(:,1);z1Org(:,2)];
uOrg2 = [z2Org(:,1);z2Org(:,2)];
uOrg3 = [z3Org(:,1);z3Org(:,2)];
uOrg4 = [z4Org(:,1);z4Org(:,2)];

uNew1 = [z1New(:,1);z1New(:,2)];
uNew2 = [z2New(:,1);z2New(:,2)];
uNew3 = [z3New(:,1);z3New(:,2)];
uNew4 = [z4New(:,1);z4New(:,2)];

[yOrg1,xOrg1] = hist(uOrg1,unique(uOrg1));
[yOrg2,xOrg2] = hist(uOrg2,unique(uOrg2));
[yOrg3,xOrg3] = hist(uOrg3,unique(uOrg3));
[yOrg4,xOrg4] = hist(uOrg4,unique(uOrg4));

[yNew1,xNew1] = hist(uNew1,unique(uNew1));
[yNew2,xNew2] = hist(uNew2,unique(uNew2));
[yNew3,xNew3] = hist(uNew3,unique(uNew3));
[yNew4,xNew4] = hist(uNew4,unique(uNew4));









