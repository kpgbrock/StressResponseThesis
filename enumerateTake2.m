function [results] = enumerateTake2()

load phycfull;
load enumerateY2HChoose3;
x = enumerateY2HChoose3;
z5 = find(phycfull.pHRange==5);
z7 = find(phycfull.pHRange==7);

results = (phycfull.charge(z5,x(:,1)).*phycfull.charge(z5,x(:,2)))./(phycfull.seqLength(x(:,1)).*phycfull.seqLength(x(:,2)))' ...
    + (phycfull.charge(z5,x(:,1)).*phycfull.charge(z5,x(:,3)))./(phycfull.seqLength(x(:,1)).*phycfull.seqLength(x(:,3)))' ...
    + (phycfull.charge(z5,x(:,3)).*phycfull.charge(z5,x(:,2)))./(phycfull.seqLength(x(:,3)).*phycfull.seqLength(x(:,2)))' ...
    - (phycfull.charge(z7,x(:,1)).*phycfull.charge(z7,x(:,2)))./(phycfull.seqLength(x(:,1)).*phycfull.seqLength(x(:,2)))' ...
    - (phycfull.charge(z7,x(:,1)).*phycfull.charge(z7,x(:,3)))./(phycfull.seqLength(x(:,1)).*phycfull.seqLength(x(:,3)))' ...
    - (phycfull.charge(z7,x(:,3)).*phycfull.charge(z7,x(:,2)))./(phycfull.seqLength(x(:,3)).*phycfull.seqLength(x(:,2)))';

results = sortrows([x results'],4);
end