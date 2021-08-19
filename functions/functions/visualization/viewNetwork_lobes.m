function viewNetwork_lobes(visual)

% setup

[mymap_conn, mymap_diff] = create_colourmap(visual);

% Make group averages and differences (all regions)

[groupadjmat,diff_groupadjmat,numGroups,numDiff,numSources,filt_freqs] = group_averages(visual);

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
    
    % generate averages for lobes
    lobe_length=length(new_lobes);
    group_lobes=zeros(lobe_length,lobe_length,numGroups);
    nan_group_nodes = group_nodes; 
    nan_group_nodes(find(group_nodes == 0)) = NaN;
    for gg=1:numGroups
        for i=1:lobe_length
            for j=1:lobe_length
                temp = nan_group_nodes(new_lobes{i},new_lobes{j},gg) ;
%                 group_lobes(i,j,gg)=mean(temp(:),'omitnan'); % average
                group_lobes(i,j,gg)=sum(temp(:),'omitnan'); % sum
            end
        end
    end
    group_lobes(isnan(group_lobes)) = 0; 

    if numDiff ~= 0
        diff_lobes=zeros(lobe_length,lobe_length,numDiff);
        nan_diff_nodes=diff_nodes; 
        nan_diff_nodes(find(diff_nodes == 0)) = NaN;
        for dd=1:numDiff
            for i=1:lobe_length
                for j=1:lobe_length
                    temp = nan_diff_nodes(new_lobes{i},new_lobes{j},dd) ;
%                     diff_lobes(i,j,dd)=mean(temp(:),'omitnan'); % average
                    diff_lobes(i,j,dd)=sum(temp(:),'omitnan'); % sum
                end
            end
        end
    end
    diff_lobes(isnan(diff_lobes)) = 0;

    % display plots
    figure('Name', [num2str(filt_freqs(visual.fb_interest,1)),'-',num2str(filt_freqs(visual.fb_interest,2)),' Hz']);
    clim = max(group_lobes(:),[],'omitnan');
    for gg=1:numGroups
    subplot(1,numGroups,gg);
    imagesc(group_lobes(:,:,gg), [0 clim]);
    title([visual.group_names{gg},' Network']);
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
            title(['Difference Network: ',visual.group_names{dd},' VS ',visual.group_names{dd+1}]);
            colorbar
            colormap(mymap_diff);
            set(gca,'XTick',[1:lobe_length],'XTickLabel',labels,...
            'YTick',[1:lobe_length],'YTickLabel',labels)
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
    
    % create a weighted mask (weighted_sig_nodes)
    temp_sig_nodes1 = nan(size(thresholded_saliences));
    temp_sig_nodes1(thresholded_saliences > 0) = thresholded_saliences(thresholded_saliences > 0);
    temp_sig_nodes1 = temp_sig_nodes1/sum(temp_sig_nodes1(:),'omitnan');
    temp_sig_nodes2 = nan(size(thresholded_saliences));
    temp_sig_nodes2(thresholded_saliences < 0) = thresholded_saliences(thresholded_saliences < 0);
    temp_sig_nodes2 = temp_sig_nodes2/sum(temp_sig_nodes2(:),'omitnan');
    weighted_sig_nodes = cat(3, temp_sig_nodes1, temp_sig_nodes2);

    network_name={'Positive','Negative'};
    
    for n=1:2
        
        for gg=1:numGroups
%             group_nodes(:,:,gg)=groupadjmat(:,:,visual.fb_interest,gg).*sig_nodes(:,:,n); % non weighted
            group_nodes(:,:,gg)=groupadjmat(:,:,visual.fb_interest,gg).*weighted_sig_nodes(:,:,n); % weighted using saliences
        end

        if numDiff~=0
            for dd=1:numDiff
%                 diff_nodes(:,:,dd)=diff_groupadjmat(:,:,visual.fb_interest,dd).*sig_nodes(:,:,n); % non weighted
                diff_nodes(:,:,dd)=diff_groupadjmat(:,:,visual.fb_interest,dd).*weighted_sig_nodes(:,:,n); % weighted using saliences
            end
        end

        % generate averages for lobes
        lobe_length=length(new_lobes);
        group_lobes=zeros(lobe_length,lobe_length,numGroups);
        nan_group_nodes=group_nodes; 
        nan_group_nodes(find(group_nodes == 0)) = NaN;
        for gg=1:numGroups
            for i=1:lobe_length
                for j=1:lobe_length
                    temp = nan_group_nodes(new_lobes{i},new_lobes{j},gg) ;
%                     group_lobes(i,j,gg)=mean(temp(:),'omitnan');     % average               
                    group_lobes(i,j,gg)=sum(temp(:),'omitnan'); % sum
                end
            end
        end
        group_lobes(isnan(group_lobes)) = 0; 

        if numDiff ~= 0
            diff_lobes=zeros(lobe_length,lobe_length,numDiff);
            nan_diff_nodes=diff_nodes; 
            nan_diff_nodes(find(diff_nodes == 0)) = NaN;
            for dd=1:numDiff
                for i=1:lobe_length
                    for j=1:lobe_length
                        temp = nan_diff_nodes(new_lobes{i},new_lobes{j},dd) ;
%                         diff_lobes(i,j,dd)=mean(temp(:),'omitnan'); % average
                        diff_lobes(i,j,dd)=sum(temp(:),'omitnan'); % sum
                    end
                end
            end
        end
        diff_lobes(isnan(diff_lobes)) = 0;

        % display plots
        figure('Name', [num2str(filt_freqs(visual.fb_interest,1)),'-',num2str(filt_freqs(visual.fb_interest,2)),' Hz']);
        clim = max(group_lobes(:),[],'omitnan');
        for gg=1:numGroups
        subplot(1,numGroups,gg);
        imagesc(group_lobes(:,:,gg), [0 clim]);
        title([visual.group_names{gg},' ',network_name{n},' Network']);
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
                title([network_name{n},' Difference Network: ',visual.group_names{dd},' VS ',visual.group_names{dd+1}]);
                colorbar
                colormap(mymap_diff);
                set(gca,'XTick',[1:lobe_length],'XTickLabel',labels,...
                'YTick',[1:lobe_length],'YTickLabel',labels)
            end
        end
        
        % number of nodes
        nnodes = nnz(any(sig_nodes));
        % number of edges
        nedges = nnz(sig_nodes)/2;
         
        fprintf(['\nPLS ',network_name{n},' Network:\n']);
        fprintf(['\t No. Nodes: ',num2str(nnodes),'\n\t No. Edges: ', num2str(nedges),'\n']);

    end
    
end

end