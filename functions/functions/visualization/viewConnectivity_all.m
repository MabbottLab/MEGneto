function viewConnectivity_all(visual)

% setup

[mymap_conn, mymap_diff] = create_colourmap(visual);

% Make group averages and differences (all regions)

[groupadjmat,diff_groupadjmat,numGroups,numDiff,~,filt_freqs] = group_averages(visual);

%  View connectivity for all regions

figure('Name', [num2str(filt_freqs(visual.fb_interest,1)),'-',num2str(filt_freqs(visual.fb_interest,2)),' Hz']);
temp = squeeze(groupadjmat(:,:,visual.fb_interest,:));
temp(temp == 1) = NaN;
clim = max(temp(:),[],'omitnan');
for gg=1:numGroups
subplot(1,numGroups,gg);
imagesc(groupadjmat(:,:,visual.fb_interest,gg),[0 clim]);
title(visual.group_names{gg});
colorbar
colormap(mymap_conn);
end

if numDiff~=0
    figure('Name', [num2str(filt_freqs(visual.fb_interest,1)),'-',num2str(filt_freqs(visual.fb_interest,2)),' Hz']);
    temp = squeeze(diff_groupadjmat(:,:,visual.fb_interest,:));
    temp(temp == 0) = NaN;
    clim = [abs(min(temp(:),[],'omitnan')) max(temp(:),[],'omitnan')] ;
    for dd=1:numDiff
    subplot(1,numDiff,dd);
    imagesc(diff_groupadjmat(:,:,visual.fb_interest,dd),[-max(clim) max(clim)]); 
    title(['Difference: ',visual.group_names{dd},' VS ',visual.group_names{dd+1}]);
    colorbar
    colormap(mymap_diff);
    end
end

end