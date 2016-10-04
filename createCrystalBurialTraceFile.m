function [dist2,pdbseq,xyz,bb] = createCrystalBurialTraceFile(pdb,seq)

% Get data from pdb file
data = pdbread(pdb);

% Vectors containing cartesian coordinates for all atoms and for the
% backbone alpha carbon atoms respectively
xyz = zeros(length(data.Model(1).Atom),3);
bb = zeros(length(seq),3);
pdbseq = [];

resCtr = 1;
for i=1:length(data.Model(1).Atom)
    
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

% Find the center of mass of the entire protein for each coordinate
xcm = mean(xyz(:,1));
ycm = mean(xyz(:,2));
zcm = mean(xyz(:,3));

% Calculate the r^2 distance of amino acids
dist2 = (bb(:,1)-xcm).^2 + (bb(:,2)-ycm).^2 + (bb(:,3)-zcm).^2;

end
