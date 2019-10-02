function viewConnectivity_lobes(visual)

% setup

[mymap_conn, mymap_diff] = create_colourmap(visual);

% Make group averages and differences (all regions)

[groupadjmat,diff_groupadjmat,numGroups,numDiff,numSources,filt_freqs] = group_averages(visual);

% View connectivity for 'Lobes' - average regions into lobes

% define regions in each lobe
frontal=[1:28,69:70];               % F = frontal
insular_cingulategyri=[29:36];      % I = insular and cingulate gyri
temporal=[37:42,55:56,79:90];       % T = temporal
occipital=[43:54];                  % O = occipital
parietal=[57:68];                   % P = parietal
central_structure=[71:78];          % S = subcortical
posterior_fossa=[91:116];           % C = cerebellar regions
if numSources==116
    lobes={frontal,insular_cingulategyri,temporal,occipital,parietal,central_structure,posterior_fossa};
    labels={'F','I','T','O','P','S','C'};
elseif numSources==90
    lobes={frontal,insular_cingulategyri,temporal,occipital,parietal,central_structure};
    labels={'F','I','T','O','P','S'};
else
    error('No lobe groupings are specified.');
end

% find where 'lobes' are located in connectiviy matrix (their index)
all_regions=sort(cell2mat(lobes));
for i=1:size(lobes,2)
    index=[];
    lobes_spread=lobes{i};
    for j=1:size(lobes_spread,2)
        index=cat(2,index,find(all_regions==lobes_spread(j)));
    end
    new_lobes{i}=index;
end

% generate averages for lobes
lobe_length=length(new_lobes);
group_lobes=zeros(lobe_length,lobe_length,numGroups);
nan_grougroupadjmat=groupadjmat; % set diagonal to NaN, so that the mean isn't biased.
nan_grougroupadjmat(find(groupadjmat == 1)) = NaN;
for gg=1:numGroups
    for i=1:lobe_length
        for j=1:lobe_length
            temp = nan_grougroupadjmat(new_lobes{i},new_lobes{j},visual.fb_interest,gg);
            group_lobes(i,j,gg)=mean(temp(:),'omitnan');
        end
    end
end

if numDiff ~= 0
    diff_lobes=zeros(lobe_length,lobe_length,numDiff);
    nan_diff_groupadjmat=diff_groupadjmat; % set diagonal to NaN, so that the mean isn't biased.
    nan_diff_groupadjmat(find(diff_groupadjmat == 0)) = NaN;
    for dd=1:numDiff
        for i=1:lobe_length
            for j=1:lobe_length
                temp = nan_diff_groupadjmat(new_lobes{i},new_lobes{j},visual.fb_interest,dd) ;
                diff_lobes(i,j,dd)=mean(temp(:),'omitnan');
            end
        end
    end
end

% display plots
figure('Name', [num2str(filt_freqs(visual.fb_interest,1)),'-',num2str(filt_freqs(visual.fb_interest,2)),' Hz']);
clim = max(group_lobes(:),[],'omitnan');
for gg=1:numGroups
subplot(1,numGroups,gg);
imagesc(group_lobes(:,:,gg), [0 clim]);
title(visual.group_names{gg});
colorbar
colormap(mymap_conn);
set(gca,'XTick',[1:lobe_length],'XTickLabel',labels,...
    'YTick',[1:lobe_length],'YTickLabel',labels)
end

if numDiff~=0
    figure('Name', [num2str(filt_freqs(visual.fb_interest,1)),'-',num2str(filt_freqs(visual.fb_interest,2)),' Hz']);
    clim = [abs(min(diff_lobes(:),[],'omitnan')) max(diff_lobes(:),[],'omitnan')] ;
    for dd=1:numDiff
        subplot(1,numDiff,dd);
        imagesc(diff_lobes(:,:,dd),[-max(clim) max(clim)]); 
        title(['Difference: ',visual.group_names{dd},' VS ',visual.group_names{dd+1}]);
        colorbar
        colormap(mymap_diff);
        set(gca,'XTick',[1:lobe_length],'XTickLabel',labels,...
        'YTick',[1:lobe_length],'YTickLabel',labels)
    end
end



end