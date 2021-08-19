function viewSaliences_all(visual)

% setup
[~, mymap_diff] = create_colourmap(visual);

% PLS data
load(visual.pls_path);

% THRESHOLD SALIENCES
thresh_percent = 0.5;
%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%

% get frequency bands
load(visual.group_paths{1},'p');
filt_freqs=p.filt_freqs;

pos_saliences=thresholded_saliences.*(thresholded_saliences > 0); % positive saliences only
neg_saliences=thresholded_saliences.*(thresholded_saliences < 0); % negative saliences only

salience = cat(3,pos_saliences,neg_saliences);

% display plots
figure('Name', [num2str(filt_freqs(visual.fb_interest,1)),'-',num2str(filt_freqs(visual.fb_interest,2)),' Hz']);
plot_titles = {'Positive Saliences','Negative Saliences'};
clim=[abs(min(neg_saliences(:))),max(pos_saliences(:))];
for ss=1:2
    
    subplot(1,2,ss);
    imagesc(salience(:,:,ss),[-1*(max(clim)) max(clim)]);
    title(plot_titles{ss});
    colorbar
    colormap(mymap_diff);
    
    
    % number of nodes
    nnodes = nnz(any(salience(:,:,ss)));
    % number of edges
    nedges = nnz(salience(:,:,ss))/2;

    fprintf(['\nPLS ',plot_titles{ss},' Network:\n']);
    fprintf(['\t No. Nodes: ',num2str(nnodes),'\n\t No. Edges: ', num2str(nedges),'\n']);
    
end



end