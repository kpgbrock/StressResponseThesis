% Make plots for self-aggregation energy, native energy, and average
% hydrophobicity for yeast proteins in a certain size range found in the
% cytosol

load yCytoExp;
energy = zeros(456,1);
selfagg = zeros(456,1);
avghydro = zeros(456,1);
ct = 1;
abundance = zeros(456,1);
for i=1:length(yCytoExp.abundance)
    fasta = yCytoExp.sequence{i};
    if (length(fasta)>=100) && (length(fasta)<=300)
        [~,~,~,emin,~,~,~,seq] = pholder(fasta);
        energy(ct) = emin;
        avghydro(ct) = mean(seq);
        [~,~,~,emin,~,~,~,~] = pholderAgg(fasta,fasta);
        selfagg(ct) = emin;
        abundance(ct) = yCytoExp.abundance(ct);
        ct = ct+1;
    end   
end

figure;
plot(abundance,avghydro,'*');
xlabel('Abundance');
ylabel('Average Hydrophobicity');
title('Abundance vs Hydrophobicity of Sequences 100-300');
print -djpeg abundancevhydro.jpg;

figure;
plot(avghydro,selfagg,'*');
xlabel('Average hydrophobicity');
ylabel('Self-Aggregation Energy');
title('Self Aggregation versus Sequence Hydrophobicity');
print -djpeg selfaggvhydro.jpg;

figure;
plot(abundance,energy,'*');
xlabel('Abundance');
ylabel('Native Minimum Energy');
title('Abundance Vs Minimized Energy');
print -djpeg abundancevenergy.jpg;

figure(energy,selfagg,'*');
xlabel('Native Minimum Energy');
ylabel('Self-aggregation energy');
title('Native Vs Self-Aggregation Energy');
print -djpeg energyvagg.jpg;
