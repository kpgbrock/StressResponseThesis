function [results] = btAndHydroEnergyCorrelation(numTimes)

orgSeq = 'MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG';
L = length(orgSeq);

hydroEnergyAA = zeros(numTimes,2);
hydroEnergyAB = zeros(numTimes,2);
hydroEnergyA = zeros(numTimes,1);
hydroEnergyB = zeros(numTimes,1);

btCorrAA = zeros(numTimes,1);
btCorrAB = zeros(numTimes,1);
btCorrOrgAAOrgAB = zeros(numTimes,4);



for i=1:numTimes
    disp(i);
    seq1 = orgSeq(randperm(length(orgSeq)));
    seq2 = orgSeq(randperm(length(orgSeq)));
    
    [copt,V,~,~,~,~,phi0,seq] = pholder(seq1);
    btOrg1 = V*copt;    
    hydroEnergyA(i) = phi0*seq*btOrg1;
    
    [copt,V,~,~,~,~,phi0,seq] = pholder(seq2);
    btOrg2 = V*copt;
    hydroEnergyB(i) = phi0*seq*btOrg2;
    
    %Self self
    [copt,V,~,~,~,~,phi0,seq] = pholderAgg(seq1,seq1);
    bt = V*copt;
    
    hydroEnergyAA(i,1) = phi0*seq(1:L)'*bt(1:L);
    hydroEnergyAA(i,2) = phi0*seq((L+1):end)'*bt((L+1):end);    
    btCorrAA(i) = corr(bt(1:L),bt((L+1):end));
    btCorrOrgAAOrgAB(i,1) = corr(btOrg1,bt(1:L));
    btCorrOrgAAOrgAB(i,2) = corr(btOrg1,bt((L+1):end));
    
    % Self/nonself
    [copt,V,~,~,~,~,phi0,seq] = pholderAgg(seq1,seq2);
    bt = V*copt;
    hydroEnergyAB(i,1) = phi0*seq(1:L)'*bt(1:L);
    hydroEnergyAB(i,2) = phi0*seq((L+1):end)'*bt((L+1):end);
    btCorrAB(i) = corr(bt(1:L),bt((L+1):end));
    btCorrOrgAAOrgAB(i,3) = corr(btOrg1,bt(1:L));
    btCorrOrgAAOrgAB(i,4) = corr(btOrg2,bt((L+1):end));
        
end

results.hydroEnergyAA = hydroEnergyAA;
results.hydroEnergyAB = hydroEnergyAB;
results.btCorrAA = btCorrAA;
results.btCorrAB = btCorrAB;
results.hydroEnergyA = hydroEnergyA;
results.hydroEnergyB = hydroEnergyB;
results.btCorrOrgAAOrgAB = btCorrOrgAAOrgAB;

end