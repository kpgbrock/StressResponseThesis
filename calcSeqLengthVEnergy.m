function [lenEnergy] = calcSeqLengthVEnergy(inputfile)

data = fastaread(inputfile);
lenEnergy = zeros(length(data),2);

for i=1:length(data)
    disp(i);
    [~,~,~,emin,~,~,~,~] = pholder(data(i,1).Sequence);
    lenEnergy(i,1) = length(data(i,1).Sequence);
    lenEnergy(i,2) = emin;
end

save lenEnergy lenEnergy
end