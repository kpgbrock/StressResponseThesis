function [dataWT,dataPerm] = checkEnergyMin()

FASTA = 'yeastCytoplasmFiltered2.txt';
rawdata = fastaread(FASTA);

dataWT = zeros(length(rawdata),1);
dataPerm = zeros(length(rawdata),1);
for i=1:length(rawdata)
    seq = rawdata(i).Sequence;
    [~,~,~,emin,~,~,~,~] = pholder(seq);
    dataWT(i) = emin;
    
    seq = seq(randperm(length(seq)));
    [~,~,~,emin,~,~,~,~] = pholder(seq);
    dataPerm(i) = emin;
end
end