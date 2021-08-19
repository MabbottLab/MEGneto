function viewSaliences_anatomy(visual,nslices)

% setup
[~, mymap_diff] = create_colourmap(visual);

% PLS data
load(visual.pls_path);

% THRESHOLD SALIENCES
thresh_percent = 0.5;
%%%%%%%%%%%%%%
% temp_thresholded_saliences = nan(size(thresholded_saliences));
% temp_thresholded_saliences(thresholded_saliences > 0) = thresholded_saliences(thresholded_saliences > 0);
% salience_min = min(temp_thresholded_saliences(:),[],'omitnan');
% salience_max = max(temp_thresholded_saliences(:),[],'omitnan');
% pos_salience_threshold = ((salience_max-salience_min)*thresh_percent) + salience_min;
% 
% temp_thresholded_saliences = nan(size(thresholded_saliences));
% temp_thresholded_saliences(thresholded_saliences < 0) = thresholded_saliences(thresholded_saliences < 0);
% salience_max = min(temp_thresholded_saliences(:),[],'omitnan');
% salience_min = max(temp_thresholded_saliences(:),[],'omitnan');
% neg_salience_threshold = salience_max - ((salience_max-salience_min)*thresh_percent);
% 
% thresholded_saliences(thresholded_saliences>neg_salience_threshold & thresholded_saliences<pos_salience_threshold) = 0;
%%%%%%%%%%%%%%%

% get frequency bands
load(visual.group_paths{1},'p');
filt_freqs=p.filt_freqs;

fullPath = which('ft_preprocessing.m');
[pathstr,~,~] = fileparts(fullPath);
% import single subject MRI (Colin 27) and AAL atlas
template= ft_read_mri([pathstr,'/template/anatomy/single_subj_T1.nii']);
atlas=ft_read_atlas([pathstr,'/template/atlas/aal/ROI_MNI_V4.nii']);
% create node strength map
anatomy=template.anatomy;
pos_atlas_map=zeros(91,109,91); 
neg_atlas_map=zeros(91,109,91);  
pos_strength=strengths_und(thresholded_saliences.*(thresholded_saliences > 0)); % strength of salience or conn at each region (sum of values)
neg_strength=strengths_und(thresholded_saliences.*(thresholded_saliences < 0)); % strength of salience or conn at each region (sum of values)

for i=1:size(thresholded_saliences,1)
    index = find(atlas.tissue == i);
    pos_atlas_map(index) = pos_strength(i);
    neg_atlas_map(index) = neg_strength(i);
end

% show montage image - DEFAULT 1:80 slices
clim=max(pos_atlas_map(:));
imoverlay(anatomy,pos_atlas_map,[-1*clim clim],[],mymap_diff,[],[],nslices,[],['Positive Saliences: ',num2str(filt_freqs(visual.fb_interest,1)),'-',num2str(filt_freqs(visual.fb_interest,2)),' Hz']); % pos
clim=abs(min(neg_atlas_map(:)));
imoverlay(anatomy,neg_atlas_map,[-1*clim clim],[],mymap_diff,[],[],nslices,[],['Negative Saliences: ',num2str(filt_freqs(visual.fb_interest,1)),'-',num2str(filt_freqs(visual.fb_interest,2)),' Hz']); % neg

plot_titles = {'Positive Saliences','Negative Saliences'};
for ss= 1:2
    if ss == 1
        nan_saliences=thresholded_saliences.*(thresholded_saliences > 0); % positive saliences only
    else
        nan_saliences=thresholded_saliences.*(thresholded_saliences < 0); % negative saliences only
    end
    % number of nodes
    nnodes = nnz(any(nan_saliences));
    % number of edges
    nedges = nnz(nan_saliences)/2;

    fprintf(['\nPLS ',plot_titles{ss},' Network:\n']);
    fprintf(['\t No. Nodes: ',num2str(nnodes),'\n\t No. Edges: ', num2str(nedges),'\n']);
end

end


