% Breaks proteins into domains of size ds and performs calculations on them

function [results] = yeastNoStapVarV2(ds)

addpath '/Users/marubi/yeastGFP/';

load yCyto30All;
yc = yCyto30All;

n = length(yc.abundance);

z = (n/2)*(n-1);

AB.avgE = cell(z,1);

permApermB.avgE = cell(z,1);

ct = 1;
for i=1:n
    disp(i);
    disp('New iteration!!!!!!!');
    
    % Get first protein sequence and two permutations
    p1 = yc.sequence{i};
    d1 = breakIntoDomains(p1,ds);
    d1(end) = [];
    
    for j=(i+1):n
        
        disp(j);
        
        % Get second sequence and a permutation
        p2 = yc.sequence{j};
        d2 = breakIntoDomains(p2,ds);
        d2(end) = [];
       
        AB.avgE{ct} = zeros(length(d1),length(d2));     
        permApermB.avgE{ct} = zeros(length(d1),length(d2));
        
        for k=1:length(d1)
            for m=1:length(d2)
                dom1 = d1{k};
                dom2 = d2{m};
                
                
                dom1perm = dom1(randperm(length(dom1)));
                dom2perm = dom2(randperm(length(dom2)));
                
                % Compute scores and energies for AB, permA-B, and permA-permB
                % respectively
                [~,~,~,emin,~,~,~,~] = pholderAgg(dom1,dom2);
                
                AB.avgE{ct}(k,m) = emin;
                
                
                [~,~,~,emin,~,~,~,~] = pholderAgg(dom1perm,dom2perm);
                
                permApermB.avgE{ct}(k,m) = emin;
               
                
                
            end
        end
        % Keep track of index when using both protein sequences
        ct = ct+1;
    end
end

results.AB = AB;
results.permApermB = permApermB;
results.seqinfo = 'Domain split but no stapling';

end


