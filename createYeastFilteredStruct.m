% Create a data structure that includes identifier, sequence, and
% expression level.  Note that you have to import yeastGFPfiltered.xlsx
% into matlab using the import wizard first, which will create 'data' and
% 'textdata' variables used below.  The resulting format is a struct that
% contains fields including the yORF, the sequence, and the abundance
% measurement.

seqdata = fastaread('orf_trans.txt');
results.id = cell(length(seqdata),1);
results.sequence = cell(length(seqdata),1);
results.abundance = nan(length(seqdata),1);
for i=1:length(seqdata)
    yorf = strtok(seqdata(i,1).Header);
    for j=1:length(textdata(:,2))
        
        if strcmp(yorf,textdata(j,2))
            results.id{i} = yorf;
            results.sequence{i} = seqdata(i,1).Sequence;
            results.abundance(i) = data(j-1,7);
            break;
        end
    end
end
