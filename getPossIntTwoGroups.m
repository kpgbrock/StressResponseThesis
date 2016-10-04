function [pairs,x] = getPossIntTwoGroups(prot1,prot2)

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
   x = length(x);
end