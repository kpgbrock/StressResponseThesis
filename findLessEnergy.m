lessenergy = [];
for i=1:length(actual)
    if ~isnan(actual(i))
        if min(perm(i,:)) < actual(i)
           lessenergy = [lessenergy; i actual(i) min(perm(i,:)) mean(perm(i,:)) std(perm(i,:))]; 
        end
    end
end