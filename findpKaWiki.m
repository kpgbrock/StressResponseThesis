function [pK,residuesPos,residuesNeg] = findpKaWiki(aa)

aa = upper(aa);
pK = NaN;

residuesPos = 'RKH';
residuesNeg = 'DEYC';

% residues = 'CDEHKNRSTUWY';
% if aa == 'C'
%     pK = 8.18;
% elseif aa == 'D'
%     pK = 3.90;
% elseif aa == 'E'
%     pK = 4.07;
% elseif aa == 'H'
%     pK = 6.04;
% elseif aa == 'K'
%     pK = 10.54;
% elseif aa == 'N'
%     pK = 5.41;
% elseif aa == 'R'
%     pK = 12.48;
% elseif aa == 'S'
%     pK = 5.68;
% elseif aa == 'T'
%     pK = 5.53;
% elseif aa == 'U'
%     pK = 5.73;
% elseif aa == 'W'
%     pK = 5.885;
% elseif aa == 'Y'
%     pK = 10.46;

nAA = length(aa);
if nAA == 1
    % Actual side chain residues
    if aa == 'D'
        pK = 3.9;
    elseif aa == 'E'
        pK = 4.3;
    elseif aa == 'R'
        pK = 12.0;
    elseif aa == 'K'
        pK = 10.5;
    elseif aa == 'H'
        pK = 6.08;
    elseif aa == 'C'
        pK = 8.28;
    elseif aa == 'Y'
        pK = 10.1;
    end
else
    % C-terminus
    if strcmp(aa,'AC')
        pK = 2.35;
    elseif strcmp(aa,'CC')
        pK = 1.92;
    elseif strcmp(aa,'DC')
        pK = 1.99;
    elseif strcmp(aa,'EC')
        pK = 2.10;
    elseif strcmp(aa,'FC')
        pK = 2.20;
    elseif strcmp(aa,'GC')
        pK = 2.35;
    elseif strcmp(aa,'HC')
        pK = 1.80;
    elseif strcmp(aa,'IC')
        pK = 2.32;
    elseif strcmp(aa,'KC')
        pK = 2.16;
    elseif strcmp(aa,'LC')
        pK = 2.33;
    elseif strcmp(aa,'MC')
        pK = 2.13;
    elseif strcmp(aa,'NC')
        pK = 2.14;
    elseif strcmp(aa,'PC')
        pK = 1.95;
    elseif strcmp(aa,'QC')
        pK = 2.17;
    elseif strcmp(aa,'RC')
        pK = 1.82;
    elseif strcmp(aa,'SC')
        pK = 2.19;
    elseif strcmp(aa,'TC')
        pK = 2.09;
    elseif strcmp(aa,'VC')
        pK = 2.39;
    elseif strcmp(aa,'WC')
        pK = 2.46;
    elseif strcmp(aa,'YC')
        pK = 2.20;
        
        % N-terminus
    elseif strcmp(aa,'AN')
        pK = 9.87;
    elseif strcmp(aa,'CN')
        pK = 10.70;
    elseif strcmp(aa,'DN')
        pK = 9.90;
    elseif strcmp(aa,'EN')
        pK = 9.47;
    elseif strcmp(aa,'FN')
        pK = 9.31;
    elseif strcmp(aa,'GN')
        pK = 9.78;
    elseif strcmp(aa,'HN')
        pK = 9.33;
    elseif strcmp(aa,'IN')
        pK = 9.76;
    elseif strcmp(aa,'KN')
        pK = 9.06;
    elseif strcmp(aa,'LN')
        pK = 9.74;
    elseif strcmp(aa,'MN')
        pK = 9.28;
    elseif strcmp(aa,'NN')
        pK = 8.72;
    elseif strcmp(aa,'PN')
        pK = 10.64;
    elseif strcmp(aa,'QN')
        pK = 9.13;
    elseif strcmp(aa,'RN')
        pK = 8.99;
    elseif strcmp(aa,'SN')
        pK = 9.21;
    elseif strcmp(aa,'TN')
        pK = 9.10;
    elseif strcmp(aa,'VN')
        pK = 9.74;
    elseif strcmp(aa,'WN')
        pK = 9.41;
    elseif strcmp(aa,'YN')
        pK = 9.21;        
    end
end
end