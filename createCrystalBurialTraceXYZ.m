function [contactmap] = createCrystalBurialTraceXYZ(pdb,seq,chainID)

% Get data from pdb file
data = getpdb(pdb);

% Vectors containing cartesian coordinates for all atoms and for the
% backbone alpha carbon atoms respectively
xyz = nan(length(data.Model(1).Atom),3);
bb = zeros(length(seq),3);
pdbseq = [];

resCtr = 1;

chainStart = NaN;
chainEnd = 0;

for i=1:length(data.Model(1).Atom)
   
    if (data.Model(1).Atom(i).chainID == chainID)
        
        if isnan(chainStart)
            chainStart = i;
        end
        
        chainEnd = i;
        
        % Extract coordinates for each atom
        xyz(i,1) = data.Model(1).Atom(i).X;
        xyz(i,2) = data.Model(1).Atom(i).Y;
        xyz(i,3) = data.Model(1).Atom(i).Z;
        
        % If it's a backbone carbon, store it specially
        if strcmp('CA',data.Model(1).Atom(i).AtomName)
            bb(resCtr,:) = xyz(i,:);
            z = aminolookup('Abbreviation',data.Model(1).Atom(i).resName);
            pdbseq = [pdbseq z(1)];
            resCtr = resCtr + 1;
        end
    end
end

xyz = xyz(chainStart:chainEnd,:);

% Find the center of mass of the entire protein for each coordinate
xcm = mean(xyz(:,1));
ycm = mean(xyz(:,2));
zcm = mean(xyz(:,3));

contactmap = zeros(length(seq));
thresh = 10;
for i=1:length(seq)
    for j=1:length(seq)
        if ((xyz(i,1)-xyz(j,1))^2 + (xyz(i,2)-xyz(j,2))^2 + (xyz(i,3)-xyz(j,3))^2) < thresh
            contactmap(i,j) = 1;
        end
    end
end

% Calculate the r^2 distance of amino acids
dist2 = (bb(:,1)-xcm).^2 + (bb(:,2)-ycm).^2 + (bb(:,3)-zcm).^2;

end
