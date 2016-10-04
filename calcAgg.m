function [AA,AB,permAA,permApermA,perm1Aperm2A,ApermB,permApermB,AandBdata,data] = calcAgg()

tic
inputfile = 'yeastCytoplasmFiltered.txt';

% Number of permutations per actual sequence
np = 10;

data = fastaread(inputfile);
%N = length(data);
N = 20;
%seqmin = 100;
%seqmax = 300;

% Variables used to store actual hydropathies for AA sequences and the
% values for random permutations of each sequence
iAndj = 0.5*(N*N - N);
AA = nan(N,1);
AB = nan(iAndj,1);
permAA = nan(N,np);
permApermA = nan(N,np);
perm1Aperm2A = nan(N,np);
ApermB = nan(iAndj,np);
permApermB = nan(iAndj,np);

%hpscores = nan(N,1);
AandBdata = nan(iAndj,4);

iter = 1;
for a=1:N
    disp(a);
    Aseq = data(a,1).Sequence;
    %hpscores(a) = computeHP(Aseq);
    
    
%     [~,~,~,emin,~,~,~,~] = pholderAgg(Aseq,Aseq);
%     AA(a) = emin;
%     
%     idx = randperm(length(Aseq));
%     [~,~,~,emin,~,~,~,~] = pholderAgg(Aseq,Aseq(idx));
%     permAA(a) = emin;
%     
%     [~,~,~,emin,~,~,~,~] = pholderAgg(Aseq(idx),Aseq(idx));
%     permApermA(a) = emin;
    
    %             idx2 = randperm(length(Aseq));
    %             [~,~,~,emin,~,~,~,~] = pholderAgg(Aseq(idx),Aseq(idx2));
    %             perm1Aperm2A(a) = emin;
    
    for b = 1:a-1
        Bseq = data(b,1).Sequence;
        %if (length(Bseq) <= seqmax) && (length(Bseq) >= seqmin) && (hpscores(b) > -1) && (hpscores(b) < 0)
        
        AandBdata(iter,:) = [a,b,length(Aseq),length(Bseq)];
        
        [~,~,~,emin,~,~,~,~] = pholderAgg(Aseq,Bseq);
        AB(iter) = emin;
        
        for j=1:np
            idx = randperm(length(Bseq));
            [~,~,~,emin,~,~,~,~] = pholderAgg(Aseq,Bseq(idx));
            ApermB(iter,j) = emin;
        end
        %
        %                     idx2 = randperm(length(Aseq));
        %                     [~,~,~,emin,~,~,~,~] = pholderAgg(Aseq(idx2),Bseq(idx));
        %                     permApermB(iter) = emin;
        
        iter = iter + 1;
        %end
    end
    %end
end
toc

save AA AA;
save AB AB;
save permApermA permApermA;
save perm1Aperm2A perm1Aperm2A;
save permApermB permApermB;
save ApermB ApermB;
save permAA permAA;
save AandBdata AandBdata;

end

% Compute hydropathy value of sequence by summing together all
% Kyte-Doolittle hydrophathies
function score = computeHP(seq)

score = 0;
for i=1:length(seq)
    score = score + KD(seq(i));
    
end
score = score/length(seq);
end