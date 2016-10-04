function [copt,V,eps,emin,r2max,msr,phi0,seq] = pholderStapleABnewK(seq1,seq2,pairs,knew)

%Start with an amino acid sequence string denoted by the variable FASTA.
%The code computes the optimal square-distance of backbone from the center
%of the protein, and plots it, returning the optimal weights for the
%burial modes, and the matrix of modes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1.5; 
% This defines the spring stiffness between monomers in
% the chain. It is chosen so that a free Gaussian coil
% has mean square inter-monomer distance of 1, so that
% other distances are measured in units of the distance
% between two adjacent alpha carbons in the backbone.

M = 10000000; 
% This weight simply ensures that the center of mass of
% the globule remains fixed at the origin.

n1 = length(seq1);
n2 = length(seq2);
n = n1+n2; 
FASTA = strcat(seq1,seq2);
% The number of residues in the amino acid sequence.

rho0 = 250/(4*pi*(4^3)/3); 
% Estimate of monomer density in proteins
% using TIM barrel as a benchmark

r2max = (3*n/rho0/4/pi)^(2/3);
% Estimate of squared max radius of globule.

phi0 = (250/n)^(2/3)/112.5; 
% The energy scale of hydropathy.
% Here, the change in globule size with
% sequence length affects the magnitude
% of the hydropathy scale used in the
% theory since the transfer energy from
% surface to core should be roughly the same
% regardless of globule size (rmax).
%%%%%%%%text%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now construct the matrix of harmonic bonds for the polymer.
bond = zeros(n,n);
for i = 1:n,
    for j = 1:n,
        bond(i,j) = (-((j + 1) == i) - ((i + 1) == j) + 2* (i == j))*k+M;
    end
end

% Leaving the ends of the polymer free to flop around
bond(1,1) = bond(1,1) - k;
bond(n,n) = bond(n,n) - k;

% Snip bonds where sequences meet
bond(n1,n1+1) = bond(n1,n1+1)+k;
bond(n1+1,n1) = bond(n1+1,n1)+k;
bond(n1,n1) = bond(n1,n1) - k;
bond(n1+1,n1+1) = bond(n1+1,n1+1)-k;

%knew = k*10;

% Add bonds where staples are present
[r,c] = size(pairs);
for i=1:r
    n1 = pairs(i,1);
    n2 = pairs(i,2);
    bond(n1,n2) = bond(n1,n2)-knew;
    bond(n2,n1) = bond(n2,n1)-knew;
    bond(n1,n1) = bond(n1,n1)+knew;
    bond(n2,n2) = bond(n2,n2)+knew;
end

% Now turn the sequence into standard Kyte-Doolittle hydropathies
seq = zeros(n,1);
for i = 1:n,
    seq(i) = KD(FASTA(i));
end

% Now construct the Hamiltonian H and diagonalize it
H = bond + phi0 * diag(seq);
[psi, epsmat] = eig(H); 
% Psi are the eigenvectors, eps the eigenvalues

eps = diag(epsmat);
% Now construct the matrix V of squared amplitudes for the eigenfunctions

V = psi.^2;
% Now compute the optimal backbone trace using linear programming.

lb = 0*ones(n,1); 
% lower bound vector preventing c values from being negative

msr = 3*r2max/5; % mean square radius fixed in ratio to r2max

options = optimset('Display','off');
[copt,emin] = linprog(eps,V,ones(n,1)*r2max,ones(1,n),msr * n,lb,[],[],options);
optstruc = V*copt;

end