function [results] = calcKD(seq)

results = zeros(1,length(seq));
for i=1:length(seq)
    
    results(i) = KD(seq(i));
end

end