function [actualhp, permhp, hpscores, data] = calcHydropathies(inputfile)
    tic
    % Number of permutations per actual sequence
    np = 20;
    data = fastaread(inputfile);
    N = length(data);
    seqmin = 100;
    seqmax = 300;
    
    % Variables used to store actual hydropathies for AA sequences and the
    % values for random permutations of each sequence
    actualhp = zeros(N,1);
    permhp = zeros(N,np);
    hpscores = zeros(N,1);
    
    for i=1:N
       h = data(i,1).Header;
       seq = data(i,1).Sequence;
       hpscores(i) = computeHP(seq);
       if (hpscores(i) < 0) && (hpscores(i) > -1)
           %if (~isempty(regexpi(h,'cytosol', 'once')) || ~isempty(regexpi(h,'cytoplasm', 'once'))) && isempty(regexpi(h,'membrane','once'))
            if isempty(regexpi(h,'membrane','once'))   
               if ((length(seq) >= seqmin) && (length(seq) <= seqmax))
                   disp(i);
                   [~,~,~,emin,~,~,~,~] = pholder(seq);
                   actualhp(i) = emin;
                   for j=1:np
                       idx = randperm(length(seq));
                       [~,~,~,emin,~,~,~,~] = pholder(seq(idx));
                       permhp(i,j) = emin;
                   end
               else
                   actualhp(i) = NaN;
               end
           else
               actualhp(i) = NaN;
           end
       else
           actualhp(i) = NaN;
       end
    end
    toc   
    save actualhp actualhp;
    save permhp permhp;
    save hpscores hpscores;
    
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