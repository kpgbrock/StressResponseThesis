function [iOn,iOff] = getOnOffIntrxns(allIntxn,phycfull)

% Need to get the total number of interactions that 'turn on' (go from same
% to opposite attraction), 'turn off' (go from opposite to same charge from
% 7 to 5), total number of unique interactions, the number of proteins that
% may decrease repulsion at ph = 5, the number of proteins that may
% increase self repulsion at ph = 5, and the protein interactions that are
% contained entirely within the group.  It also returns the 
% protein interactions that are contained in the group.  It takes in the
% 2-vector matrix of known interaction binary pairs (in index form) and the
%list of protein indices.


% Get charges for interacting pairs at 5 and 7
z5 = find(phycfull.pHRange == 5);
z7 = find(phycfull.pHRange == 7);

% allIntxn = [0 0];
% 
% for i=1:length(protList)
%    allIntxn = [allIntxn; knownInt(knownInt(:,1) == protList(i),:)];
%    allIntxn = [allIntxn; knownInt(knownInt(:,2) == protList(i),:)];       
% end
% 
% allIntxn(1,:) = [];
allIntxn(:,3:6) = [phycfull.charge(z5,allIntxn(:,1))' phycfull.charge(z5,allIntxn(:,2))' phycfull.charge(z7,allIntxn(:,1))' phycfull.charge(z7,allIntxn(:,2))'];

% % Get unique rows of interacting pairs and store duplicates (these are
% % protein pairs that interact with each other within the group
% [allIntxn,a,~] = unique(allIntxn,'rows');
% ixDupRows = setdiff(1:size(allIntxn,1), a);
% iSameAll = allIntxn(ixDupRows,:);
% [r,~] = size(allIntxn);
% iTotal = r;


% Count how many are 'on' and how many are 'off' interacting
indOn = find((allIntxn(:,3).*allIntxn(:,4) < 0) & (allIntxn(:,5).*allIntxn(:,6) > 0));
indOff = find((allIntxn(:,3).*allIntxn(:,4) > 0) & (allIntxn(:,5).*allIntxn(:,6) < 0));

iOn = allIntxn(indOn,:);
iOff = allIntxn(indOff,:);

end