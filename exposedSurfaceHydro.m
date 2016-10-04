function [expKDresults] = exposedSurfaceHydro()

load yCytoExp;
n = length(yCytoExp.abundance);
block = 150;

expKDresults = yCytoExp;
expKDresults.avgexp = zeros(n,1);
for i=1:n
    seq = yCytoExp.sequence{i};
    
    % Break into domains
    z = ceil(length(seq)/block);
    
    KDsum = 0;
    for j=0:(z-1)
        if j == (z-1)
            fastaseq = seq((block*j+1):end);
        else
            fastaseq = seq((block*j+1):(block*(j+1)));
        end
        [copt,V,~,~,~,msr,~,hydro] = pholder(fastaseq);
        r = V*copt;
        
        % Want to find stretches on the surface (defined as having
        % residues greater than the mean radius away from the center of
        % the protein
        locs = r>(msr*1.5);
        KDsum = KDsum + sum(hydro(locs));
        
    end
    
    expKDresults.avgexp(i) = KDsum/length(seq);
    disp(i);
end
    
end