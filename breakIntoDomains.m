% Function that breaks a long sequence into domains of length up to dsize.
% It returns a cell array of the domains.

function [domains] = breakIntoDomains(seq,dsize)

n = ceil(length(seq)/dsize);
domains = cell(n,1);
for i=0:(n-2)
    domains{i+1} = seq((i*dsize+1):(i*dsize+dsize));
end
domains{n} = seq(((n-1)*dsize+1):end);
end