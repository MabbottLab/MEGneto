function viewSaliences_lobes(visual)

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

% get dimensions and frequency bands
load(visual.group_paths{1},'p');
m = matfile(p.paths.connmat(1,visual.connmetric));
[numSources, ~, ~, ~] = size(m.adjmat); % must be the same for each subject
filt_freqs=p.filt_freqs;

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

% generate averages for lobes (POSITIVE SALIENCES)
lobe_length=length(new_lobes);
pos_salience_lobes=zeros(lobe_length,lobe_length);
nan_saliences=thresholded_saliences.*(thresholded_saliences > 0); % positive saliences only
nan_saliences(find(nan_saliences == 0)) = NaN; % set zeros to NaN 
for i=1:lobe_length
    for j=1:lobe_length
        temp = nan_saliences(new_lobes{i},new_lobes{j}) ;
%         pos_salience_lobes(i,j)=mean(temp(:),'omitnan'); % average
        pos_salience_lobes(i,j)=sum(temp(:),'omitnan'); % sum
    end
end
pos_salience_lobes(isnan(pos_salience_lobes)) = 0;


% generate averages for lobes (NEGATIVE SALIENCES)
lobe_length=length(new_lobes);
neg_salience_lobes=zeros(lobe_length,lobe_length);
nan_saliences=thresholded_saliences.*(thresholded_saliences < 0); % negative saliences only
nan_saliences(find(nan_saliences == 0)) = NaN; % set zeros to NaN 
for i=1:lobe_length
    for j=1:lobe_length
        temp = nan_saliences(new_lobes{i},new_lobes{j}) ;
%         neg_salience_lobes(i,j)=mean(temp(:),'omitnan'); % average
        neg_salience_lobes(i,j)=sum(temp(:),'omitnan'); % sum
    end
end
neg_salience_lobes(isnan(neg_salience_lobes)) = 0;

salience_lobes = cat(3,pos_salience_lobes,neg_salience_lobes);

% display plots
figure('Name', [num2str(filt_freqs(visual.fb_interest,1)),'-',num2str(filt_freqs(visual.fb_interest,2)),' Hz']);
plot_titles = {'Positive Saliences','Negative Saliences'};
clim=[abs(min(neg_salience_lobes(:))),max(pos_salience_lobes(:))];
for ss=1:2
    subplot(1,2,ss);
    imagesc(salience_lobes(:,:,ss),[-1*(max(clim)) max(clim)]);
    title(plot_titles{ss});
    colorbar
    colormap(mymap_diff);
    set(gca,'XTick',[1:lobe_length],'XTickLabel',labels,...
        'YTick',[1:lobe_length],'YTickLabel',labels)
    
    
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