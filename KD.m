function [ kd ] = KD( AA )
%Takes amino acid letter code as its input and returns Kyte-Doolittle
%hydropathy.

kd = 0;

if AA =='I', 
    kd = 4.5;
end

if AA =='X', 
    kd = -0.4;
end

if AA == 'V',
    kd = 4.2;
end

if AA == 'L',
    kd = 3.8;
end

if AA == 'F',
    kd = 2.8;
end

if AA == 'C',
    kd = 2.5;
end

if AA == 'M',
    kd = 1.9;
end

if AA == 'A',
    kd = 1.8;
end

if AA == 'G',
    kd = -0.4;
end

if AA == 'T',
    kd = -0.7;
end

if AA == 'W',
    kd = -0.9;
end

if AA == 'S',
    kd = -0.8;
end

if AA == 'Y',
    kd = -1.3;
end

if AA == 'P',
    kd = -1.6;
end

if AA == 'H',
    kd = -3.2;
end

if AA == 'E',
    kd = -3.5;
end

if AA == 'Q',
    kd = -3.5;
end

if AA == 'D',
    kd = -3.5;
end

if AA == 'N',
    kd = -3.5;
end

if AA == 'K',
    kd = -3.9;
end

if AA == 'R',
    kd = -4.5;
end


end

