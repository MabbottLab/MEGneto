function Step4c_DynamicFC(config, pid, visitnum)
% DynamicFC calculates specified FC metric using sliding window analysis
% and ft_connectivityanalysis
%
% Inputs:
%       connmethod: string, see ft_connectivityanalysis (ex: wpli_debiased)
%       bandavgmethod: string (ex: max)
%       toi: time of interest in seconds, double (ex: [0 1])
%       winsize: size of sliding window in seconds, float  
%       stepsize: how much time window slides over, in seconds, float
%       freqbands: array of doubles containing frequency ranges
%
% Returns: 
%       .mat containing connectivity matrices (window number x ROI x ROI x
%       frequency band)
%       .png plotting metric of ROI x ROI over time windows for each frequency 
%
% Notes:
%       - uses all trials/does not separate by conditions
%       - if win_size/step_size aren't a perfect multiple of the full time
%       window, truncates at last multiple

%% SETUP
    if exist('visitnum', 'var')
        this_output = [config.meta.project_path '/' config.meta.analysis_name '/' pid '/' sprintf('ses-%.2d', visitnum)]; % indicate subject-specific output folder path
    else
        this_output = [config.meta.project_path '/' config.meta.analysis_name '/' pid]; 
    end

    load([this_output '/out_struct.mat'])
    load([this_output '/step3_data_roi.mat'])

%% RUN DFC ANALYSIS
    connDFC        = [];
    connDFC.dimord = 'chan_chan_freq';
    connDFC.label  = data_roi.label;
    
    t_start   = config.step4c.toi(1); 
    t_end     = config.step4c.toi(2);
    winsize  = config.step4c.winsize;
    stepsize = config.step4c.stepsize;
    
    if stepsize > winsize
        warning('Chosen step size is larger than window size');
        
    end 
    
    %%% FOR EACH FREQUENCY BAND --------------------------------------------
    for fq = 1:size(config.step4c.freqbands,1) 
        % filter data
        cfg           = [];
        cfg.bpfilter  = 'yes';
        cfg.bpfreq    = config.step4c.freqbands(fq,:);
        data_roi_filt = ft_preprocessing(cfg, data_roi);
        
        % variables that get reused 
        win_num = 1; 
        f_low   = config.step4c.freqbands(fq,1);
        f_high  = config.step4c.freqbands(fq,2);
        
        % 1/f_low = time for one cycle of the wave in seconds
        if 1/f_low > winsize | 1/f_high > winsize
            warning(sprintf('Chosen window size captures less than one cycle of frequency band %d - %d Hz',f_low,f_high))
        end 
        
        fprintf('Frequency band: %d - %d Hz \n',f_low,f_high);
        
        %%% FOR EACH TIME WINDOW -------------------------------------------
        for t = t_start:stepsize:(t_end-winsize)
            % select data based on smaller window 
            cfg               = [];
            cfg.latency       = [t t+ winsize];
            data_roi_filt_toi = ft_selectdata(cfg,data_roi_filt); 
            
            % keep track of centre of window 
            connDFC.win_centre{fq}(win_num) = t + winsize/2; 
            
            % CALCULATE CONNECTIVITY 
            fprintf('Connectivity calculations for window: %.3f to %.3f s \n',t,t + winsize)
            cfg             = []; % set up config for connectivity calculation
            cfg.method      = 'mtmfft';
            cfg.output      = 'powandcsd';
            cfg.channel     = 'all';
            cfg.trials      = 'all';
            cfg.keeptrials  = 'yes';
            cfg.taper       = 'hanning';
            cfg.foilim      = config.step4c.freqbands(fq,:);
            freq_filt       = ft_freqanalysis(cfg, data_roi_filt_toi); 

            %%% COMPUTE CONNECTIVITY METRIC -------------------------------
            cfg             = []; % set up config for computing connectivity metric
            cfg.method      = config.step4c.connmethod;
            this_conn       = ft_connectivityanalysis(cfg, freq_filt); % compute connectivity metric

            % RESHAPE INTO SOURCE X SOURCE CONN MAT AND STORE -------------
            this_conn        = ft_checkdata(this_conn, 'cmbstyle', 'full'); % the one wpli debiased is a matrix nsensor x nsensor x nfreq matrix 
            connDFC.freq(fq) = mean(this_conn.freq);
          
            % collapse frequencies within band
            if strcmp(config.step4c.bandavgmethod, 'max') % default to max across band
                connDFC.(sprintf('%sspctrm',cfg.method))(win_num,:,:,fq) = ...
                    squeeze(max(this_conn.(sprintf('%sspctrm', cfg.method)),[], 3));
            else % otherwise, use mean
                connDFC.(sprintf('%sspctrm',cfg.method))(win_num,:,:,fq) = ...
                    squeeze(mean(this_conn.(sprintf('%sspctrm', cfg.method)), 3));
            end            
            win_num = win_num + 1; 
        end 
        
        % -----------------------------------------------------------------------------
        % SAVING ROI X ROI MATRICES OVER TIME AS FIGURE FOR CURRENT
        % FREQUENCY BAND
        % -----------------------------------------------------------------------------
        matrix   = connDFC.(sprintf("%sspctrm",config.step4c.connmethod));
        matrix_fq = matrix(:,:,:,fq); % selecting data for this freq band
        n_win    = size(matrix, 1);
            
        % ROI labels
        labels = this_conn.label;
        connDFC.label = labels; %re-assigning b/c ft_checkdata re-orders ROIs
        
        % window time labels for plots
        t_init = t_start;
        t_fin  = t_init + winsize; 
        
        % store max and min color values for uniform colorbar
        colorLims = [];
        
        % prevent figures from popping up 
        set(0, 'DefaultFigureVisible','off'); 
        
        % move onto rows of 4 cols if more than 4 time windows
        if n_win > 4
            if mod(n_win,4)==0 % if n_win is divisible by 4
                n_row = n_win/4;
            else
                n_row = floor(n_win/4) + 1;
            end
            n_col = 4; 
        else
            n_row = 1;
            n_col = n_win;
        end 
        
        % create tiled layout figure for heatmaps
        fig = figure(fq); 
        tcl = tiledlayout(fig,n_row,n_col,'TileSpacing','Compact'); 
        h   = gobjects(n_row,n_col); 
        
        % plot heatmap for each window 
        for w = 1:n_win            
            nexttile(tcl);
            
            % heatmap at window "w"
            h(w) = heatmap(labels,labels,squeeze(matrix_fq(w,:,:)));
            
            %%% formatting 
            if mod(w+3,4)~=0 % if the index+3 is not divisible by 4 (i.e. not in the first col)
                h(w).YDisplayLabels = nan(size(h(w).YDisplayLabels)); % remove ylabels
            end
            
            colorbar off
            h(w).NodeChildren(3).TickLabelInterpreter = 'none'; % prevent letter subscripts for ROI names
            set(h(w),'GridVisible','off','FontSize',6); 
            colorLims=[colorLims; h(w).ColorLimits];
            
            % label matrix with time window 
            title_text = compose([sprintf("%.3f to %.3f s",t_init,t_fin)]); 
            h(w).Title = '\fontsize{10}' + title_text;
            t_init     = t_init + stepsize;
            t_fin      = t_init + winsize;
            
        end 
        
        title(tcl,sprintf('Dynamic Functional Connectivity For Frequency %d-%d Hz',f_low,f_high));
        
        % make colorbar constant for figure 
        globalColorLim = [min(colorLims(:,1)), max(colorLims(:,2))];
        set(h(1), 'ColorLimits', globalColorLim)
        ax = axes(tcl,'visible','off','Colormap',h(1).Colormap,'CLim',globalColorLim);
        cb             = colorbar(ax);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
        cb.Layout.Tile = 'East';

        % figure labeled with freq band 
        saveas(fig,[this_output sprintf('dfc_freq_%d_%d.png',f_low,f_high)])
               
    end 
    
    save([this_output '/step4c_connDFC.mat'], 'connDFC'); 
    
end