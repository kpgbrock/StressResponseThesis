function [results] = getOutsideGO()

load yCyt;
load phycfull;
load godata;

N = length(yCyt.id);
ngo = length(godata.pIdx);
results.goVCyto = zeros(ngo,3);
results.goVCytoOnInd = cell(ngo,1);
results.goVCytoOffInd = cell(ngo,1);
goPairs = nchoosek(1:ngo,2);
results.goVgoGroupInd = goPairs;
results.goVgoNumSame = zeros(length(goPairs),1);

% Get charges for interacting pairs at 5 and 7
z5 = find(phycfull.pHRange == 5);
z7 = find(phycfull.pHRange == 7);

%How does protein group interact with all other proteins in cytosol, 
%compared to within itself?
for i=1:ngo
   proteins = godata.pIdx{i};
   refIndx = (1:N)';
   refIndx(proteins) = [];
   [p,q] = meshgrid(proteins, refIndx);
   pairs = [p(:) q(:)];
   z = [phycfull.charge(z5,pairs(:,1))' phycfull.charge(z5,pairs(:,2))' phycfull.charge(z7,pairs(:,1))' phycfull.charge(z7,pairs(:,2))'];
   results.goVCytoOnInd{i} = pairs((z(:,1).*z(:,2) < 0) & (z(:,3).*z(:,4) > 0),:);
   results.goVCytoOffInd{i} = pairs((z(:,1).*z(:,2) > 0) & (z(:,3).*z(:,4) < 0),:);
   [r1,~] = size(results.goVCytoOnInd{i});
   [r2,~] = size(results.goVCytoOffInd{i});
   results.goVCyto(i,:) = [r1 r2 length(pairs)];

end

disp('Out of first loop');
%Interaction of starvation stress with punctae proteins

%Remove soluble (2nd) set from 1st punctae set and rerun calculations on just this set

%For each 2 GO groups:  look at all possible interactions, P1 in 1 and P2 
%in 2, and see which pairs are the most interactive

for i=1:length(goPairs)
   prot1 = godata.pIdx{goPairs(i,1)};
   prot2 = godata.pIdx{goPairs(i,2)};

   [p,q] = meshgrid(prot1,prot2);
   pairs = [p(:) q(:)];
   
   % Take out equivalent pairings, like 1 2 and 2 1
   [r,~] = size(pairs);
   for j=1:r
      if (pairs(j,1) > pairs(j,2))
          pairs(j,:) = [pairs(j,2) pairs(j,1)];
      end
   end
   pairs = unique(pairs,'rows');
   
   % Store number of proteins in common between the two groups, and the
   % total number of interactions when the self-self interactions have been
   % removed   
   x = find(pairs(:,1)==pairs(:,2));
   pairs(x,:) = [];
   
   [r,~] = size(pairs);
   results.goVgoNumSame(i) = length(x);  
%    if i==1
%        disp(pairs);
%    end
   z = [phycfull.charge(z5,pairs(:,1))' phycfull.charge(z5,pairs(:,2))' phycfull.charge(z7,pairs(:,1))' phycfull.charge(z7,pairs(:,2))'];
   zOn = pairs((z(:,1).*z(:,2) < 0) & (z(:,3).*z(:,4) > 0),:);
   zOff = pairs((z(:,1).*z(:,2) > 0) & (z(:,3).*z(:,4) < 0),:);
   [r1,~] = size(zOn);
   [r2,~] = size(zOff);
   results.goVgo(i,:) = [r1 r2 r];
end

end