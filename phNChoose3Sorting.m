function [z] = phNChoose3Sorting()

load enumerateY2HChoose3;
load phycfull;
N = length(enumerateY2HChoose3);
z = zeros(N,3);

for i=1:N
    if mod(i,100000)==0
        disp(i);
    end
    protIntrxns = nchoosek(enumerateY2HChoose3(i,:),2);
    data = isingModelIntrxns(protIntrxns, 7,5);
    z(i,1) = data.eTotalNeutral;
    z(i,2) = data.eTotalAcidic;
    z(i,3) = data.dE;
    
end

z = [enumerateY2HChoose3 z];
z = sortrows(z,6);
end