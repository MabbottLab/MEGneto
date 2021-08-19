function interpret_PLS(result, numSubj, numCond, numGroups, numTrials, numSources, numSamples, numDataPoints, network, fb_interest)

%%
if  network.pls.method == 1 || network.pls.method == 4
    if network.pls.method == 4
        [task_range,~]=size(result.TBv{1,1});
        u=result.u(:,1:task_range);
        v=result.TBv{1,1};
        s=result.s(1:task_range);
        boot_ratio=(result.boot_result.compare_u(:,1:task_range));
        probability=result.perm_result.sprob(1:task_range);
        brain_scores=result.usc; % FIX - make sure the correct ranges are selected
    else
        u=result.u;
        v=result.v;
        s=result.s;
        boot_ratio=result.boot_result.compare_u;
        probability=result.perm_result.sprob;
        brain_scores=result.usc;
    end
    clear temp_contrast contrast pval
    for lv=1:network.lv_num
        if strcmp(network.plsmetric,'lcmv') % time series
            saliences = nan(numSources, numSamples, numTrials);
            bootstrapratio = nan(numSources, numSamples, numTrials);
            for tt = 1:numTrials
                col = (((tt-1)*numDataPoints))+(1:numDataPoints);  

                % for saliences
                tmp_matrix = reshape(u(col,lv),numSamples,numSources);
                saliences(:,:,tt) = tmp_matrix';

                % for bootstrap ratio
                tmp_matrix = reshape(boot_ratio(col,lv),numSamples,numSources);
                bootstrapratio(:,:,tt) = tmp_matrix';

                % threshold salince by bootstrap ratio
                bootratio_threshold(tt,1:2) = quantile(boot_ratio(col,lv),[network.q_tail (1-network.q_tail)]);
                thrpos = bootratio_threshold(tt,2)
                thrneg = bootratio_threshold(tt,1)
                neg_thresholded_saliences(:,:,tt) = saliences(:,:,tt).*(bootstrapratio(:,:,tt) < thrneg);
                pos_thresholded_saliences(:,:,tt) = saliences(:,:,tt).*(bootstrapratio(:,:,tt) > thrpos);
                thresholded_saliences = neg_thresholded_saliences + pos_thresholded_saliences; % no overlapping
            end
        else % connectivity
            saliences = nan(numSources, numSources);
            bootstrapratio = nan(numSources, numSources);
            col = 1:numDataPoints;  

            % for saliences
%             disp(size(col));
%             disp(size(lv));
            tmp_triangle = squareform(u(col,lv), 'tomatrix');
            tmp_matrix = tmp_triangle;
%             size(tmp_matrix)
%             size(saliences)
            saliences(:,:) = tmp_matrix;

            % for bootstrap ratio
            tmp_triangle = squareform(boot_ratio(col,lv), 'tomatrix');    
            tmp_matrix = tmp_triangle;    
            bootstrapratio(:,:) = tmp_matrix;
            
            % threshold salience by bootstrap ratio
            bootratio_threshold(1:2) = quantile(boot_ratio(col,lv),[network.q_tail (1-network.q_tail)]);
            thrpos = bootratio_threshold(2);
            thrneg = bootratio_threshold(1);
            neg_thresholded_saliences(:,:) = saliences(:,:).*(bootstrapratio(:,:) < thrneg);
            pos_thresholded_saliences(:,:) = saliences(:,:).*(bootstrapratio(:,:) > thrpos);
            thresholded_saliences = neg_thresholded_saliences + pos_thresholded_saliences; % no overlapping
            
        end

        % other info
        temp_contrast = v(:,lv); % re-arrange contrast
        for gg = 1:numGroups
            contrast(gg,:) = temp_contrast((1:numCond) + numCond*(gg-1));
        end 
        pval = probability(lv); % significance of latent variable

        
      % VISUALIZATION
      
        % contrast bar graph
        if strcmp(network.plsmetric,'lcmv') == 0 % time series
            hf = figure('Name', [num2str(network.filt_freqs(fb_interest,1)),'-',num2str(network.filt_freqs(fb_interest,2)),' Hz']);
        else
            hf = figure;
        end
        bar(contrast)
        grid on;
        if numCond > 1 && numGroups > 1
            legend(network.cond_names,1)
        end
        title(['LV ',num2str(lv),' Contrast [p-value = ', num2str(pval),' ]']);
        if numGroups == 1 && numCond > 1 % just use conditions names instead
            set(gca,'XtickLabel',network.cond_names,'Box','off','GridLineStyle','none');
        else
            set(gca,'XtickLabel',network.group_names,'Box','off','GridLineStyle','none');
        end
        print(hf, [network.file_out,'/PLS_',num2str(network.filt_freqs(fb_interest,1)),'-',num2str(network.filt_freqs(fb_interest,2)),'Hz_lv',num2str(lv),'_contrast.png'], '-dpng', '-r600');

        % bootstrap ratios histogram
        if strcmp(network.plsmetric,'lcmv') == 0 % time series
            hf = figure('Name', [num2str(network.filt_freqs(fb_interest,1)),'-',num2str(network.filt_freqs(fb_interest,2)),' Hz']);
        else
            hf = figure;
        end
        [tmp1,tmp2] = hist(bootstrapratio(:,:),500);
        h00 = bar(tmp2,sum(tmp1,2),0.7,'stack','FaceColor',[0.4 0.4 0.4]);
        hold on;
        bin_select = find(tmp2<thrneg);
        h02 = bar(tmp2(bin_select),sum(tmp1(bin_select,:),2),0.7,'stack','FaceColor',[155 175 228]/256,'EdgeColor','b');
        hold on;
        bin_select = find(tmp2>thrpos);
        h01 = bar(tmp2(bin_select),sum(tmp1(bin_select,:),2),0.7,'stack','FaceColor',[237 139 135]/256,'EdgeColor','r');
        vline(thrpos,'r--');
        vline(thrneg,'b--');
        set(gca,'Box','off');
        title(['LV ',num2str(lv),' BSR Threshold: ',num2str(thrneg),' & ',num2str(thrpos)]);
        xlabel('Bootstrap Ratios (BSR)'); %,'FontSize',14,'FontWeight','bold');
        ylabel('Counts'); %,'FontSize',14,'FontWeight','bold');        
        print(hf, [network.file_out,'/PLS_',num2str(network.filt_freqs(fb_interest,1)),'-',num2str(network.filt_freqs(fb_interest,2)),'Hz_lv',num2str(lv),'_BSRhist.png'], '-dpng', '-r600');
        
        if strcmp(network.plsmetric,'lcmv') % time series
         % boot strap ratio time series (per source) 
            Source = 100;
            highlighted_brs=zeros(numSources,numSamples);
            highlighted_brs(bootstrapratio(:,:,1) < thrneg)=1;
            highlighted_brs(bootstrapratio(:,:,1) > thrpos)=1;
            highlighted_brs=highlighted_brs.*bootstrapratio(:,:,1);
            figure;
            plot(linspace(0,numSamples/srate,numSamples),bootstrapratio(Source,:,1),'Color',[0.4 0.4 0.4]);
            hold on;
            plot(linspace(0,numSamples/srate,numSamples),highlighted_brs(Source,:,1),'g');
            l1=refline(0,thrpos);
            l2=refline(0,thrneg);
            set(l1,'Color','r','LineStyle','--');
            set(l2,'Color','b','LineStyle','--');
            set(gca,'Box','off');
            xlabel('Time (s)');
            ylabel('Bootstrap Ratios')
            title (['LV ',num2str(lv),' Bootstrap Ratio Timeseries']);
        % group average/correlations with bootstrap threshold time series (per source?)
        
        
            
        % bootstrap thresholded salience image over time
            fullPath = which('ft_preprocessing.m');
            [pathstr,~,~] = fileparts(fullPath);
            % import single subject MRI (Colin 27) and AAL atlas
            template= ft_read_mri([pathstr,'/template/anatomy/single_subj_T1.nii']);
            atlas=ft_read_atlas([pathstr,'/template/atlas/aal/ROI_MNI_V4.nii']);
            % create node strength map
            anatomy=template.anatomy;
            atlas_map=atlas.tissue;  
            % show montage image - DEFAULT 1:80 slices
            mymap= brewermap(100,'*RdBu');
            imoverlay_timeseries(anatomy,atlas_map,thresholded_saliences,mymap); % FIX
            
        else % connectivity
        % bootstrap thresholded salience image (salience strength)
            fullPath = which('ft_preprocessing.m');
            [pathstr,~,~] = fileparts(fullPath);
            % import single subject MRI (Colin 27) and AAL atlas
            template= ft_read_mri([pathstr,'/template/anatomy/single_subj_T1.nii']);
            atlas=ft_read_atlas([pathstr,'/template/atlas/aal/ROI_MNI_V4.nii']);
            % create node strength map
            anatomy=template.anatomy;
            pos_atlas_map=zeros(91,109,91); 
            neg_atlas_map=zeros(91,109,91);  
            pos_salience_strength=strengths_und(pos_thresholded_saliences(:,:)); % strength of salience at each region (sum of saliences)
            neg_salience_strength=strengths_und(neg_thresholded_saliences(:,:)); % strength of salience at each region (sum of saliences)
            for i=1:numSources 
                index = find(atlas.tissue == i);
                pos_atlas_map(index) = pos_salience_strength(i);
                neg_atlas_map(index) = neg_salience_strength(i);
            end
            atlas_map = pos_atlas_map + neg_atlas_map; % sum of saliences for each region (pos and neg network superimposed)
            
            % show montage image - DEFAULT 1:80 slices
            mymap= brewermap(100,'*RdBu');
            if nnz(atlas_map) == 0
                disp('Salience was insignificant');
%                 clim=[-1,1];
            else
                clim=[abs(min(atlas_map(:))),abs(max(atlas_map(:)))]                         
                [~,~,hf] = imoverlay(anatomy,atlas_map,[-1*(max(clim)) max(clim)],[],mymap);         
                print(hf, [network.file_out,'/PLS_',num2str(network.filt_freqs(fb_interest,1)),'-',num2str(network.filt_freqs(fb_interest,2)),'Hz_lv',num2str(lv),'_BSRthresholded_saliences.png'], '-dpng', '-r600'); % FIX printing image with wrong dimensions

            end 

        end
        
        % brain scores for each subject in each group (shows how each inidivual expresses the LV)
        if network.crossover == 1
            temp_BS = zeros(max(numSubj),numGroups,numCond);
            temp_numSubj=[0 numSubj];
            for gg= 1:numGroups
                for cc = 1:numCond
                    BS_grouprange = (temp_numSubj(gg+1)*(cc-1))+(sum(temp_numSubj(1:gg))*numCond)+(1:temp_numSubj(gg+1));
                    temp_BS(1:numSubj(gg),gg,cc)=brain_scores(BS_grouprange,lv);
                end
            end
            hf = figure('Position',[200 200 1000 500]);
            for cc = 1:numCond
                subplot(numCond,1,cc)
                bar(temp_BS(:,:,cc));
                title(['LV ',num2str(lv),' Brain Scores: ',network.cond_names{cc}]);
                legend(network.group_names);
            end
            xlabel('Subjects in Group');
        else
            temp_BS = zeros(max(numSubj),numCond,numGroups);
            temp_numSubj=[0 numSubj];
            for gg= 1:numGroups
                for cc = 1:numCond
                    BS_grouprange = (temp_numSubj(gg+1)*(cc-1))+(sum(temp_numSubj(1:gg))*numCond)+(1:temp_numSubj(gg+1));
                    temp_BS(1:numSubj(gg),cc,gg)=brain_scores(BS_grouprange,lv);
                end
            end
            hf = figure('Position',[200 200 1000 500]);
            for gg = 1:numGroups
                subplot(numGroups,1,gg)
                bar(temp_BS(:,:,gg));
                title(['LV ',num2str(lv),' Brain Scores: ',network.group_names{gg}]);
                if numCond > 1
                    legend(network.cond_names);
                end
            end
            xlabel('Subjects in Group');
        end
        print(hf, [network.file_out,'/PLS_',num2str(network.filt_freqs(fb_interest,1)),'-',num2str(network.filt_freqs(fb_interest,2)),'Hz_lv',num2str(lv),'_brainscores.png'], '-dpng', '-r600');

        % save data
        save([network.file_out,'/PLS_',num2str(network.filt_freqs(fb_interest,1)),'-',num2str(network.filt_freqs(fb_interest,2)),'Hz_lv',num2str(lv),'.mat'],'thresholded_saliences','-mat');

    end
    
    % latent variable (LV) effect size (singular values)
    hf = figure;
    bar(s);
    ylim([0 1]);
    title('Latent Variable Singular Values');
    ylabel('Effect size of each LV (%)');
    xlabel('Latent Variables (LV)');
    text((1:length(s))-0.26,double(s)+0.03,num2str(s));
    print(hf, [network.file_out,'/PLS_',num2str(network.filt_freqs(fb_interest,1)),'-',num2str(network.filt_freqs(fb_interest,2)),'Hz_LVeffect.png'], '-dpng', '-r600');
        
    % latent variable (LV) p-value 
    hf = figure;
    bar(probability);
    ylim([0 1]);
    title('Latent Variable P-values');
    ylabel('P-value of each LV (%)');
    xlabel('Latent Variables (LV)');
    text((1:length(probability))-0.2,double(probability)+0.03,num2str(probability));
    print(hf, [network.file_out,'/PLS_',num2str(network.filt_freqs(fb_interest,1)),'-',num2str(network.filt_freqs(fb_interest,2)),'Hz_LVpval.png'], '-dpng', '-r600');

end

% --------------------------------------------------------------------------
    
if  network.pls.method == 3 || network.pls.method == 4
    if network.pls.method == 4
        [behav_range,~]=size(result.TBv{1,2});
        u=result.u(:,task_range+1:task_range+behav_range);
        v=result.TBv{1,2};
        s=result.s(task_range+1:task_range+behav_range);
        boot_ratio=(result.boot_result.compare_u(:,task_range+1:task_range+behav_range));
        probability=result.perm_result.sprob(task_range+1:task_range+behav_range);
        brain_scores=result.usc; % FIX - make sure the correct ranges are selected
    else
        u=result.u;
        v=result.v;
        s=result.s;
        boot_ratio=result.boot_result.compare_u;
        probability=result.perm_result.sprob;
        brain_scores=result.usc;
    end
    clear temp_contrast contrast pval
    numBehave = size(network.pls.stacked_behavdata,2);
    for lv=1:network.lv_num
        if strcmp(network.plsmetric,'lcmv') % time series
            saliences = nan(numSources, numSamples, numTrials, numBehave);
            bootstrapratio = nan(numSources, numSamples, numTrials, numBehave);
            for bb=1:numBehave
                for tt = 1:numTrials
                    col = (((tt-1)*numDataPoints))+(1:numDataPoints);  
                    behav=((lv*numBehave)-(numBehave-1))+(bb-1); 

                    % for saliences
                    tmp_matrix = reshape(u(col,behav),numSamples,numSources);
                    saliences(:,:,tt,bb) = tmp_matrix';

                    % for bootstrap ratio
                    tmp_matrix = reshape(boot_ratio(col,behav),numSamples,numSources);
                    bootstrapratio(:,:,tt,bb) = tmp_matrix';

                    % threshold salince by bootstrap ratio
                    
                    bootratio_threshold(tt,1:2,bb) = quantile(boot_ratio(col,behav),[network.q_tail (1-network.q_tail)]);

                    neg_thresholded_saliences(:,:,tt,bb) = saliences(:,:,tt,bb).*(bootstrapratio(:,:,tt,bb) < bootratio_threshold(tt,1,bb));
                    pos_thresholded_saliences(:,:,tt,bb) = saliences(:,:,tt,bb).*(bootstrapratio(:,:,tt,bb) > bootratio_threshold(tt,2,bb));
                    thresholded_saliences = neg_thresholded_saliences + pos_thresholded_saliences; % no overlapping
                end            
                % other info
                pval(bb) = probability(behav); % significance of latent variable
                temp_contrast = v(:,lv); % contrast
                % reorganize correlation matrix
                for gg = 1:numGroups
                    contrast(bb,gg) = temp_contrast((2*gg-1)+(bb-1)); % correlation
                end 
            end            
        else % connectivity
            saliences = nan(numSources, numSources, numBehave);
            bootstrapratio = nan(numSources, numSources, numBehave);
            contrast_title = ['LV ',num2str(lv),' Correlation  [p-value = ']; % title used further down (needed to be out of loop)
            for bb=1:numBehave
                col = 1:numDataPoints;  
                behav=((lv*numBehave)-(numBehave-1))+(bb-1);

                % for saliences
                tmp_triangle = squareform(u(col,behav), 'tomatrix');
                tmp_matrix = tmp_triangle;
                saliences(:,:,bb) = tmp_matrix;

                % for bootstrap ratio
                tmp_triangle = squareform(boot_ratio(col,behav), 'tomatrix');    
                tmp_matrix = tmp_triangle;    
                bootstrapratio(:,:,bb) = tmp_matrix;
          
                % threshold salience by bootstrap ratio
                bootratio_threshold(1:2,bb) = quantile(boot_ratio(col,behav),[network.q_tail (1-network.q_tail)]);

                neg_thresholded_saliences(:,:,bb) = saliences(:,:,bb).*(bootstrapratio(:,:,bb) < bootratio_threshold(1,bb));
                pos_thresholded_saliences(:,:,bb) = saliences(:,:,bb).*(bootstrapratio(:,:,bb) > bootratio_threshold(2,bb));
                thresholded_saliences = neg_thresholded_saliences + pos_thresholded_saliences; % no overlapping

                % other info
                pval(bb) = probability(behav); % significance of latent variable
                temp_contrast = v(:,lv); % correlation
                LVeffect(:,bb) = s(behav:numBehave:length(s));
                LVpval(:,bb) = probability(behav:numBehave:length(probability));
                
                % reorganize correlation matrix
                if numGroups == 1 && numCond > 1 
                    for cc = 1:numCond
                        contrast(bb,cc) = temp_contrast(((cc*numBehave)-(numBehave-1))+(bb-1)); 
                    end 
                    % make title for contrast matrix
                    contrast_title = [contrast_title, num2str(pval(bb)),' (',network.pls.behavdata{bb},') '];
                elseif numGroups > 1 && numCond == 1
                    for gg = 1:numGroups
                        contrast(bb,gg) = temp_contrast(((gg*numBehave)-(numBehave-1))+(bb-1)); 
                    end 
                     % make title for contrast matrix
                    contrast_title = [contrast_title,num2str(pval(bb)),' (',network.pls.behavdata{bb},') '];
                elseif numGroups > 1 && numCond > 1
                    for gg = 1:numGroups
                        for cc = 1:numCond
                            contrast(bb,cc,gg) = temp_contrast(((gg*(numBehave*numCond))-((numBehave*numCond)-1))+(cc-1)+(bb-1)); 
                        end 
                         % make title for contrast matrix % FIX
                        contrast_title = [contrast_title,num2str(pval(bb)),' (',network.pls.behavdata{bb},') '];
                    end
                end
            end
        end
      
      % VISUALIZATION
          
        % correlation bar graph
        if strcmp(network.plsmetric,'lcmv') == 0 % time series
            figure('Name', [num2str(network.filt_freqs(fb_interest,1)),'-',num2str(network.filt_freqs(fb_interest,2)),' Hz']);
        else
            figure;
        end
        if numCond > 1 && numGroups == 1
            bar(contrast);   
            legend(network.cond_names,1);
            grid on
            title([contrast_title,']']); 
            set(gca,'XtickLabel',network.pls.behavdata, 'Box','off','GridLineStyle','none');
        elseif numCond == 1 && numGroups > 1
            bar(contrast);
            legend(network.group_names,1);
            grid on
            title([contrast_title,']']); 
            set(gca,'XtickLabel',network.pls.behavdata, 'Box','off','GridLineStyle','none');
        elseif numCond > 1 && numGroups > 1
            for gg = 1:numGroups
                subplot(1,numGroups,gg);
                bar(contrast(:,:,gg));
                legend(network.cond_names,1);
                grid on
                title(network.group_names{gg});
                set(gca,'XtickLabel',network.pls.behavdata, 'Box','off','GridLineStyle','none');
                % text([],[],contrast_title)                % FIX - add p-value to image
            end
        end
        print(hf, [network.file_out,'/PLS_',num2str(network.filt_freqs(fb_interest,1)),'-',num2str(network.filt_freqs(fb_interest,2)),'Hz_lv',num2str(lv),'_contrast.png'], '-dpng', '-r600');
        

        % bootstrap ratios histogram
        if strcmp(network.plsmetric,'lcmv') == 0 % time series
            figure('Name', [num2str(network.filt_freqs(fb_interest,1)),'-',num2str(network.filt_freqs(fb_interest,2)),' Hz']);
        else
            figure;
        end
        for bb=1:numBehave
            subplot(1,numBehave,bb);
            thrneg = bootratio_threshold(1,bb);
            thrpos = bootratio_threshold(2,bb);
            [tmp1,tmp2] = hist(bootstrapratio(:,:,bb),500);
            h00 = bar(tmp2,sum(tmp1,2),0.7,'stack','FaceColor',[0.4 0.4 0.4]);
            hold on;
            bin_select = find(tmp2<thrneg);
            h02 = bar(tmp2(bin_select),sum(tmp1(bin_select,:),2),0.7,'stack','FaceColor',[155 175 228]/256,'EdgeColor','b');
            hold on;
            bin_select = find(tmp2>thrpos);
            h01 = bar(tmp2(bin_select),sum(tmp1(bin_select,:),2),0.7,'stack','FaceColor',[237 139 135]/256,'EdgeColor','r');
            vline(thrpos,'r--');
            vline(thrneg,'b--');
            set(gca,'Box','off');
            title([char(network.pls.behavdata(bb)),' BSR Threshold: ',num2str(thrneg),' & ',num2str(thrpos)]);
            xlabel('Bootstrap Ratios'); %,'FontSize',14,'FontWeight','bold');
            ylabel('Counts'); %,'FontSize',14,'FontWeight','bold'); 
            hold on;
        end
        print(hf, [network.file_out,'/PLS_',num2str(network.filt_freqs(fb_interest,1)),'-',num2str(network.filt_freqs(fb_interest,2)),'Hz_lv',num2str(lv),'_BSRhist.png'], '-dpng', '-r600');
         
        if strcmp(network.plsmetric,'lcmv') % time series
         % boot strap ratio time series (per source) 
            Source = 100;
            highlighted_brs=zeros(numSources,numSamples);
            highlighted_brs(bootstrapratio(:,:,1) < thrneg)=1;
            highlighted_brs(bootstrapratio(:,:,1) > thrpos)=1;
            highlighted_brs=highlighted_brs.*bootstrapratio(:,:,1);
            figure;
            plot(linspace(0,numSamples/srate,numSamples),bootstrapratio(Source,:,1),'Color',[0.4 0.4 0.4]);
            hold on;
            plot(linspace(0,numSamples/srate,numSamples),highlighted_brs(Source,:,1),'g');
            l1=refline(0,thrpos);
            l2=refline(0,thrneg);
            set(l1,'Color','r','LineStyle','--');
            set(l2,'Color','b','LineStyle','--');
            set(gca,'Box','off');
            xlabel('Time (s)');
            ylabel('Bootstrap Ratios')
            title ('Bootstrap Ratio Timeseries');
        % group average/correlations with bootstrap threshold time series (per source?)
            
        % bootstrap thresholded salience image over time
            fullPath = which('ft_preprocessing.m');
            [pathstr,~,~] = fileparts(fullPath);
            % import single subject MRI (Colin 27) and AAL atlas
            template= ft_read_mri([pathstr,'/template/anatomy/single_subj_T1.nii']);
            atlas=ft_read_atlas([pathstr,'/template/atlas/aal/ROI_MNI_V4.nii']);
            % create node strength map
            anatomy=template.anatomy;
            atlas_map=atlas.tissue;  
            % show montage image - DEFAULT 1:80 slices
            mymap= brewermap(100,'*RdBu');
            imoverlay_timeseries(anatomy,atlas_map,thresholded_saliences,mymap); % FIX
            
        else % connectivity
        % bootstrap thresholded salience image
            fullPath = which('ft_preprocessing.m');
            [pathstr,~,~] = fileparts(fullPath);
            % import single subject MRI (Colin 27) and AAL atlas
            template= ft_read_mri([pathstr,'/template/anatomy/single_subj_T1.nii']);
            atlas=ft_read_atlas([pathstr,'/template/atlas/aal/ROI_MNI_V4.nii']);
            % create node strength map
            anatomy=template.anatomy;
            pos_atlas_map=zeros(91,109,91); 
            neg_atlas_map=zeros(91,109,91); 
            for bb=1:numBehave
                pos_salience_strength=strengths_und(pos_thresholded_saliences(:,:,bb)); % strength of salience at each region (sum of saliences)
                neg_salience_strength=strengths_und(neg_thresholded_saliences(:,:,bb)); % strength of salience at each region (sum of saliences)
                for i=1:numSources 
                    index = find(atlas.tissue == i);
                    pos_atlas_map(index) = pos_salience_strength(i);
                    neg_atlas_map(index) = neg_salience_strength(i);
                end
                atlas_map = pos_atlas_map + neg_atlas_map; % sum of saliences for each region (pos and neg network superimposed)
                % show montage image - DEFAULT 1:80 slices
                mymap= brewermap(100,'*RdBu');
                clim=[abs(min(atlas_map(:))),abs(max(atlas_map(:)))];
                [hf,~,~] = imoverlay(anatomy,atlas_map,[-1*(max(clim)) max(clim)],[],mymap,[],[],[],[],[char(network.pls.behavdata(bb)),': BSR Thresholded Salience Strength']);
                print(hf, [network.file_out,'/PLS_',num2str(network.filt_freqs(fb_interest,1)),'-',num2str(network.filt_freqs(fb_interest,2)),'Hz_lv',num2str(lv),'_BSRthresholded_saliences_',char(network.pls.behavdata(bb)),'.png'], '-dpng', '-r600'); % FIX printing image with wrong dimensions
            end           
        end
        
        % brain scores for each subject in each group (shows how each inidivual expresses the LV)
        if network.crossover == 1
            temp_BS = zeros(max(numSubj),numGroups,numCond);
            temp_numSubj=[0 numSubj];
            for gg= 1:numGroups
                for cc = 1:numCond
                    BS_grouprange = (temp_numSubj(gg+1)*(cc-1))+(sum(temp_numSubj(1:gg))*numCond)+(1:temp_numSubj(gg+1));
                    temp_BS(1:numSubj(gg),gg,cc)=brain_scores(BS_grouprange,lv);
                end
            end
            hf = figure('Position',[200 200 1000 500]);
            for cc = 1:numCond
                subplot(numCond,1,cc)
                bar(temp_BS(:,:,cc));
                title(['LV ',num2str(lv),' Brain Scores: ',network.cond_names{cc}]);
                legend(network.group_names);
            end
            xlabel('Subjects in Group');
        else
            temp_BS = zeros(max(numSubj),numCond,numGroups);
            temp_numSubj=[0 numSubj];
            for gg= 1:numGroups
                for cc = 1:numCond
                    BS_grouprange = (temp_numSubj(gg+1)*(cc-1))+(sum(temp_numSubj(1:gg))*numCond)+(1:temp_numSubj(gg+1));
                    temp_BS(1:numSubj(gg),cc,gg)=brain_scores(BS_grouprange,lv);
                end
            end
            hf = figure('Position',[200 200 1000 500]);
            for gg = 1:numGroups
                subplot(numGroups,1,gg)
                bar(temp_BS(:,:,gg));
                title(['LV ',num2str(lv),' Brain Scores: ',network.group_names{gg}]);
                if numCond > 1
                    legend(network.cond_names);
                end
            end
            xlabel('Subjects in Group');
        end
        print(hf, [network.file_out,'/PLS_',num2str(network.filt_freqs(fb_interest,1)),'-',num2str(network.filt_freqs(fb_interest,2)),'Hz_lv',num2str(lv),'_brainscores.png'], '-dpng', '-r600');
        
        % save data
        save([network.file_out,'/PLS_',num2str(network.filt_freqs(fb_interest,1)),'-',num2str(network.filt_freqs(fb_interest,2)),'Hz_lv',num2str(lv),'.mat'],'thresholded_saliences','-mat');

    end  
    
        
    % latent variable (LV) effect size (singular values)
    hf = figure;
    bar(LVeffect);
    title('Latent Variable Singular Values');
    ylabel('Effect size of each LV');
    xlabel('Latent Variables (LV)');
    legend(network.pls.behavdata);
%     text((1:length(s))-0.26,double(s)+0.03,num2str(s));
    print(hf, [network.file_out,'/PLS_',num2str(network.filt_freqs(fb_interest,1)),'-',num2str(network.filt_freqs(fb_interest,2)),'Hz_LVeffect.png'], '-dpng', '-r600');

    % latent variable (LV) p-value 
    hf = figure;
    bar(LVpval);
    ylim([0 1]);
    title('Latent Variable P-values');
    ylabel('P-value of each LV (%)');
    xlabel('Latent Variables (LV)');
    legend(network.pls.behavdata);
%     text((1:length(probability))-0.2,double(probability)+0.03,num2str(probability));
    print(hf, [network.file_out,'/PLS_',num2str(network.filt_freqs(fb_interest,1)),'-',num2str(network.filt_freqs(fb_interest,2)),'Hz_LVpval.png'], '-dpng', '-r600');

    
        
end


% close all;

end