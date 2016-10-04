function [results] = mhmcGroupY2H(nProts,maxIters)

    load y2hdata;
    load yCyt;
    load phycfull;
    
    N = length(y2hdata.idxIndProtNS);
    %protlist = 1:N;
    protGroup = randperm(N,nProts);
    x = isingModelIntrxns(nchoosek(protGroup,2),7,5,phycfull);
    eOrg = x.dE;
    
    results.dE = zeros(maxIters,1);
    
    results.protGroup = [];
    
    for i=1:maxIters
       
        % Choose one random protein to randomly swap
        randProtID = randperm(nProts,1);
        newProtGroup = protGroup;
        protlist = 1:N;
        protlist = protlist(~ismember(protlist,protGroup));
        newProtGroup(randProtID) = protlist(randperm(length(protlist),1));
        x2 = isingModelIntrxns(nchoosek(newProtGroup,2),7,5,phycfull);
        eNew = x2.dE;
        %disp(exp(-1*(eNew-eOrg)/abs(eNew)))
        if rand(1) < exp(-1*(eNew-eOrg)*100)
            eOrg = eNew;
            protGroup = newProtGroup;
            %disp(eNew-eOrg)
        end
        
        results.dE(i) = eOrg;
    end
    results.protGroup = sort(protGroup);    
end