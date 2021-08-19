function [labels] = regions_extractinfo(regions)

%define regions of interest that can distinguish lobe group
regions_interest={7,29,85,52,60,78,101};
labels_interest={'F','I','T','O','P','S','C'};
labels={};

% differentiate between lobe groups and make labels
for i=1:length(regions)
    for j=1:7
        if find(regions{i}==regions_interest{j})~=0
            labels=cat(2,labels,labels_interest(j));
            break
        end
    end
end


return
