function viewNetwork_all(visual)

% setup

[mymap_conn, mymap_diff] = create_colourmap(visual);

% Make group averages and differences (all regions)

[groupadjmat,diff_groupadjmat,numGroups,numDiff,~,filt_freqs] = group_averages(visual);

if strcmp(visual.netmetric,'nbs')
    
    % View NBS significant nodes highlighted
    load(visual.nbs_path);

    % make binary mask for nodes in network
    sig_nodes=full(nbs.NBS.con_mat{1,1}+nbs.NBS.con_mat{1,1}'); % binary mask (significant nodes)

    for gg=1:numGroups
        group_nodes(:,:,gg)=groupadjmat(:,:,visual.fb_interest,gg).*sig_nodes;
    end

    if numDiff~=0
        for dd=1:numDiff
            diff_nodes(:,:,dd)=diff_groupadjmat(:,:,visual.fb_interest,dd).*sig_nodes;
        end
    end

    figure('Name', [num2str(filt_freqs(visual.fb_interest,1)),'-',num2str(filt_freqs(visual.fb_interest,2)),' Hz']);
    temp = group_nodes;
    temp(temp == 0 ) = NaN;
    clim = max(temp(:),[],'omitnan');
    for gg=1:numGroups
      subplot(1,numGroups,gg);
      imagesc(group_nodes(:,:,gg),[0 clim]);
      title([visual.group_names{gg},' Network']);
      colorbar
      colormap(mymap_conn);
    end

    if numDiff~=0
        figure('Name', [num2str(filt_freqs(visual.fb_interest,1)),'-',num2str(filt_freqs(visual.fb_interest,2)),' Hz']);
        temp = diff_nodes;
        temp(temp == 0 ) = NaN;
        clim = [abs(min(temp(:),[],'omitnan')) max(temp(:),[],'omitnan')];
        for dd=1:numDiff
            subplot(1,numDiff,dd);
            imagesc(diff_nodes(:,:,dd),[-max(clim) max(clim)]); 
            title(['Difference Network: ',visual.group_names{dd},' VS ',visual.group_names{dd+1}]);
            colorbar
            colormap(mymap_diff);
        end
    end

    % number of nodes
    nnodes = nnz(any(sig_nodes));
    % number of edges
    nedges = nnz(sig_nodes)/2;

    fprintf(['\nNBS Network:\n']);
    fprintf(['\t No. Nodes: ',num2str(nnodes),'\n\t No. Edges: ', num2str(nedges),'\n']);

elseif strcmp(visual.netmetric,'pls')
    
    % View PLS significant nodes highlighted in Connectivity Matrices
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

    % create binary mask (sig_nodes)
    sig_nodes = cat(3,thresholded_saliences > 0,thresholded_saliences < 0);
    network_name={'Positive','Negative'};
    
    for i=1:2
       
        for gg=1:numGroups
            group_nodes(:,:,gg)=groupadjmat(:,:,visual.fb_interest,gg).*sig_nodes(:,:,i);
        end

        if numDiff~=0
            for dd=1:numDiff
                diff_nodes(:,:,dd)=diff_groupadjmat(:,:,visual.fb_interest,dd).*sig_nodes(:,:,i);
            end
        end

        figure('Name', [num2str(filt_freqs(visual.fb_interest,1)),'-',num2str(filt_freqs(visual.fb_interest,2)),' Hz']);
        temp = group_nodes;
        temp(temp == 0 ) = NaN;
        clim = max(temp(:),[],'omitnan');
        for gg=1:numGroups
          subplot(1,numGroups,gg);
          imagesc(group_nodes(:,:,gg),[0 clim]);
          title([visual.group_names{gg},' ',network_name{i},' Network']);
          colorbar
          colormap(mymap_conn);
        end

        if numDiff~=0
            figure('Name', [num2str(filt_freqs(visual.fb_interest,1)),'-',num2str(filt_freqs(visual.fb_interest,2)),' Hz']);
            temp = diff_nodes;
            temp(temp == 0 ) = NaN;
            clim = [abs(min(temp(:),[],'omitnan')) max(temp(:),[],'omitnan')];
            for dd=1:numDiff
            subplot(1,numDiff,dd);
            imagesc(diff_nodes(:,:,dd),[-max(clim) max(clim)]); 
            title([network_name{i},' Difference Network: ',visual.group_names{dd},' VS ',visual.group_names{dd+1}]);
            colorbar
            colormap(mymap_diff);
            end
        end
        
        % number of nodes
        nnodes = nnz(any(sig_nodes));
        % number of edges
        nedges = nnz(sig_nodes)/2;
         
        fprintf(['\nPLS ',network_name{i},' Network:\n']);
        fprintf(['\t No. Nodes: ',num2str(nnodes),'\n\t No. Edges: ', num2str(nedges),'\n']);

    end
    
end

end