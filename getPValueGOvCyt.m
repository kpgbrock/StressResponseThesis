function [results] = getPValueGOvCyt(nIts)

load yCyt;
load phycfull;
load godata;

% Get number of proteins per group
nP = zeros(length(godata.pIdx),1);
for i=1:length(godata.pIdx)
    nP(i) = length(godata.pIdx{i});
end

nPorg = nP;
nP = unique(nP);
results.nP = nP;
results.dist = zeros(length(nP),nIts);
results.percentile = zeros(length(nP),1);

% Get charges for interacting pairs at 5 and 7
z5 = find(phycfull.pHRange == 5);
z7 = find(phycfull.pHRange == 7);

for i=1:length(nP)
    disp(i);
    for j=1:nIts
        z = (1:length(yCyt.id))';
        z = z(randperm(length(z)));
        
        zP = z(1:nP(i));
        zCyt = z((nP(i)+1):end);
        
        [p,q] = meshgrid(zP, zCyt);
        pairs = [p(:) q(:)];
        
        [r,~] = size(pairs);
        z = [phycfull.charge(z5,pairs(:,1))' phycfull.charge(z5,pairs(:,2))' phycfull.charge(z7,pairs(:,1))' phycfull.charge(z7,pairs(:,2))'];
        [rOn,~] = size(pairs((z(:,1).*z(:,2) < 0) & (z(:,3).*z(:,4) > 0),:));
        %[rOff,~] = size(pairs((z(:,1).*z(:,2) > 0) & (z(:,3).*z(:,4) < 0),:));
        
        results.dist(i,j) = rOn/r;        
    end    
end

gocyto = getOutsideGO();
actFrac = gocyto.goVCyto(:,1)./gocyto.goVCyto(:,3);
for i=1:length(actFrac)
    dist1 = results.dist(nP == length(godata.pIdx{i}),:);
    results.percentile(i) = length(find(dist1 < actFrac(i)))/nIts;
end
end