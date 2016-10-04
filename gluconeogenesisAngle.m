function [results] = gluconeogenesisAngle(nTrials)

close all

load yCyt;
load phycfull;

% Negative regulation of gluconeogenesis
glucListNeg = [78; 93; 163; 428; 635; 821; 846; 1310];
glucListPos = [402; 596];

rGlucNeg = isingModelIntrxns(nchoosek(glucListNeg,2), 7,5, phycfull);
rGlucPos = isingModelIntrxns(nchoosek(glucListPos,2),7,5,phycfull);

results.actNeutralNeg = rGlucNeg.eTotalNeutral;
results.actAcidicNeg = rGlucNeg.eTotalAcidic;
results.actDENeg = rGlucNeg.dE;
results.actFracNegNeg = length(find((rGlucNeg.eComponentsAcidic - rGlucNeg.eComponentsNeutral) < 0))/length(rGlucNeg.eComponentsAcidic);
results.eTotalNeutralNeg = zeros(nTrials,1);
results.eTotalAcidicNeg = zeros(nTrials,1);
results.dENeg = zeros(nTrials,1);
results.fracNegNeg = zeros(nTrials,1);

results.actNeutralPos = rGlucPos.eTotalNeutral;
results.actAcidicPos = rGlucPos.eTotalAcidic;
results.actDEPos = rGlucPos.dE;
results.actFracNegPos = length(find((rGlucPos.eComponentsAcidic - rGlucPos.eComponentsNeutral) < 0))/length(rGlucPos.eComponentsAcidic);
results.eTotalNeutralPos = zeros(nTrials,1);
results.eTotalAcidicPos = zeros(nTrials,1);
results.dEPos = zeros(nTrials,1);
results.fracNegPos = zeros(nTrials,1);

results.rpETotalNeutralNeg = zeros(nTrials,1);
results.rpETotalAcidicNeg = zeros(nTrials,1);
results.rpDENeg = zeros(nTrials,1);

results.rpETotalNeutralPos = zeros(nTrials,1);
results.rpETotalAcidicPos = zeros(nTrials,1);
results.rpDEPos = zeros(nTrials,1);

% Random permutations
seqStringNeg = strjoin(yCyt.sequence(glucListNeg));
seqStringNeg(seqStringNeg == '*') = [];
seqStringNeg(seqStringNeg == ' ') = [];
seqLengthsNeg = phycfull.seqLength(glucListNeg);

seqStringPos = strjoin(yCyt.sequence(glucListPos));
seqStringPos(seqStringPos == '*') = [];
seqStringPos(seqStringPos == ' ') = [];
seqLengthsPos = phycfull.seqLength(glucListPos);

for i=1:nTrials
    NNeg = length(glucListNeg);
    protlistNeg = randperm(1804,NNeg);
    r = isingModelIntrxns(nchoosek(protlistNeg,2),7,5,phycfull);
    results.eTotalNeutralNeg(i) = r.eTotalNeutral;
    results.eTotalAcidicNeg(i) = r.eTotalAcidic;
    results.dENeg(i) = r.dE;
    results.fracNegNeg(i) = length(find((r.eComponentsAcidic - r.eComponentsNeutral) < 0))/length(r.eComponentsAcidic);
    
    %
    seqStringN = seqStringNeg(randperm(length(seqStringNeg)));
    tempdata.seq = mat2cell(seqStringN,1,seqLengthsNeg);
    tempdata.seqLength = seqLengthsNeg;
    tempdata.pHRange = phycfull.pHRange;
    tempdata.charge = zeros(length(phycfull.pHRange),length(tempdata.seq));
    
    rpNeg = getRandPermsGluco(tempdata);
    results.rpETotalNeutralNeg(i) = rpNeg.eTotalNeutral;
    results.rpETotalAcidicNeg(i) = rpNeg.eTotalAcidic;
    results.rpDENeg(i) = rpNeg.dE;
    
    seqStringP = seqStringNeg(randperm(length(seqStringPos)));
    tempdata.seq = mat2cell(seqStringP,1,seqLengthsPos);
    tempdata.seqLength = seqLengthsPos;
    tempdata.pHRange = phycfull.pHRange;
    tempdata.charge = zeros(length(phycfull.pHRange),length(tempdata.seq));
    
    rpPos = getRandPermsGluco(tempdata);
    results.rpETotalNeutralPos(i) = rpPos.eTotalNeutral;
    results.rpETotalAcidicPos(i) = rpPos.eTotalAcidic;
    results.rpDEPos(i) = rpPos.dE;
    
    %
    
    NPos = length(glucListPos);
    protlistPos = randperm(1804,NPos);
    r = isingModelIntrxns(nchoosek(protlistPos,2),7,5,phycfull);
    results.eTotalNeutralPos(i) = r.eTotalNeutral;
    results.eTotalAcidicPos(i) = r.eTotalAcidic;
    results.dEPos(i) = r.dE;
    results.fracNegPos(i) = length(find((r.eComponentsAcidic - r.eComponentsNeutral) < 0))/length(r.eComponentsAcidic);
    
end

% negNeutral = length(find(results.eTotalNeutralNeg < results.actNeutralNeg))
% negAcidic = length(find(results.eTotalAcidicNeg < results.actAcidicNeg))
% negDE = length(find(results.dENeg < results.actDENeg))
% negNumAtt = length(find(results.fracNegNeg < results.actFracNegNeg))
% 
% posNeutral = length(find(results.eTotalNeutralPos < results.actNeutralPos))
% posAcidic = length(find(results.eTotalAcidicPos < results.actAcidicPos))
% posDE = length(find(results.dEPos < results.actDEPos))
% posNumAtt = length(find(results.fracNegPos < results.actFracNegPos))

negNeutral = length(find(results.rpETotalNeutralNeg < results.actNeutralNeg))
negAcidic = length(find(results.rpETotalAcidicNeg < results.actAcidicNeg))
negDE = length(find(results.rpDENeg < results.actDENeg))

posNeutral = length(find(results.rpETotalNeutralPos < results.actNeutralPos))
posAcidic = length(find(results.rpETotalAcidicPos < results.actAcidicPos))
posDE = length(find(results.rpDEPos < results.actDEPos))


x1 = linspace(-0.06,0.14,20);
x2 = linspace(-0.01,0.01,20);
lw = 6;
fs = 24;

hold on;
subplot(3,1,1);
hist(results.eTotalNeutralNeg,x1);
axis([min(x1) max(x1) 0 5000]);
%title('Negative Regulation');
ylimits = get(gca,'YLim');
hold on;
plot([results.actNeutralNeg results.actNeutralNeg],ylimits,'k');
xlabel('\Sigma E_{Neutral}');
ylabel('Counts');
set(gca,'LineWidth',lw,'FontSize',fs);
box off;

subplot(3,1,2);
hist(results.eTotalAcidicNeg,x1);
axis([min(x1) max(x1) 0 5000]);
ylimits = get(gca,'YLim');
hold on;
plot([results.actAcidicNeg results.actAcidicNeg],ylimits,'k');
xlabel('\Sigma E_{Acidic}');
ylabel('Counts');
set(gca,'LineWidth',lw,'FontSize',fs);
box off;

subplot(3,1,3);
hist(results.dENeg,x1);
axis([min(x1) max(x1) 0 5000]);
hold on;
plot([results.actDENeg results.actDENeg],get(gca,'YLim'),'k');
xlabel('\Sigma \Delta E');
ylabel('Counts');
set(gca,'LineWidth',lw,'FontSize',fs);
box off;

figure;

subplot(3,1,1);
hold on;
hist(results.eTotalNeutralPos,x2);
axis([min(x2) max(x2) 0 5000]);
plot([results.actNeutralPos results.actNeutralPos],get(gca,'YLim'),'k');
%title('Positive Regulation');
xlabel('\Sigma E_{Neutral}');
ylabel('Counts');
set(gca,'LineWidth',lw,'FontSize',fs);
box off;


subplot(3,1,2);
hist(results.eTotalAcidicPos,x2);
axis([min(x2) max(x2) 0 5000]);
hold on;
plot([results.actAcidicPos results.actAcidicPos],get(gca,'YLim'),'k');
xlabel('\Sigma E_{Acidic}');
ylabel('Counts');
set(gca,'LineWidth',lw,'FontSize',fs);
box off;

subplot(3,1,3);
hist(results.dEPos,x2);
axis([min(x2) max(x2) 0 5000]);
hold on;
plot([results.actDEPos results.actDEPos],get(gca,'YLim'),'k');
xlabel('\Sigma \Delta E');
ylabel('Counts');
set(gca,'LineWidth',lw,'FontSize',fs);
box off;

figure;
subplot(3,1,1);
hist(results.rpETotalNeutralNeg);
hold on;
plot([results.actNeutralNeg results.actNeutralNeg],get(gca,'YLim'),'k');
axis([-0.07 0.07 0 5000])
xlabel('\Sigma E_{Neutral}');
subplot(3,1,2);
hist(results.rpETotalAcidicNeg);
hold on;
plot([results.actAcidicNeg results.actAcidicNeg],get(gca,'YLim'),'k');
axis([-0.07 0.07 0 5000])
xlabel('\Sigma E_{Acidic}');
subplot(3,1,3);
hist(results.rpDENeg);
hold on;
plot([results.actDENeg results.actDENeg],get(gca,'YLim'),'k');
axis([-0.07 0.07 0 5000])
xlabel('\Sigma \Delta E');

figure;
subplot(3,1,1);
hist(results.rpETotalNeutralPos);
hold on;
plot([results.actNeutralPos results.actNeutralPos],get(gca,'YLim'),'k');
axis([-0.004 0.002 0 5000]);
xlabel('\Sigma E_{Neutral}');
subplot(3,1,2);
hist(results.rpETotalAcidicPos);
hold on;
plot([results.actAcidicPos results.actAcidicPos],get(gca,'YLim'),'k');
axis([-0.004 0.002 0 5000]);
xlabel('\Sigma E_{Acidic}');
subplot(3,1,3);
hist(results.rpDEPos);
hold on;
plot([results.actDEPos results.actDEPos],get(gca,'YLim'),'k');
axis([-0.004 0.002 0 5000]);

end

function [data2] = getRandPermsGluco(tempdata)

% Get charge array:  findPICurveWiki
for k=1:length(tempdata.seq)
    [fcharge,~] = findPICurveWiki(tempdata.seq{k},tempdata.pHRange);
    tempdata.charge(:,k) = fcharge;
end
%

randIntrxns = nchoosek(1:length(tempdata.seq),2);


data2 = isingModelIntrxns(randIntrxns,7,5, tempdata);
%             results.histAcidic(j,i) = data2.eTotalAcidic;
%             results.histNeutral(j,i) = data2.eTotalNeutral;
%             results.histDE(j,i) = data2.dE;
end