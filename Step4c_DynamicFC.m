function Step4c_DynamicFC(config, pid, visitnum)
%
% Step4c_DynamicFC calculates specified functional connectivity metric 
% across a series of sliding windows, thereby returning dynamic functional
% connectivity rather than static (i.e., 1 value for connectivity between
% pairwise ROIs). 
%
% INPUTS:
%    > config: 
%       struct, configured in your "main" script with all analysis
%       parameters and paths. specifically, need the following fields:
%           > config.step4c.connmethod: string, see ft_connectivityanalysis (ex: wpli_debiased)
%           > config.step4c.toi: time of interest in seconds, double (ex: [0 1])
%           > config.step4c.winsize: size of sliding window in seconds, float  
%           > config.step4c.stepsize: how much time window slides over, in seconds, float
%           > config.step4c.freqbands: array of doubles containing frequency ranges
%
% Additional configuration options:
%       > config.step4c.bandavgmethod: string (ex: max)
%
% OUTPUTS: 
%       > step4c_connDFC.mat:
%           struct containing connectivity matrices (window number x ROI x ROI x
%           frequency band)
%
% Notes:
%       - uses all trials/does not separate by conditions
%       - if win_size/step_size aren't a perfect multiple of the full time
%       window, cuts off at the last window
%       ex: if toi is 0-1 s and both win_size and step_size are 0.3 s
%           windows will be: 0-0.3 s, 0.3-0.6 s, 0.6-0.9 s, leaving out
%           the last 0.1 s
%
% Last updated by: Bianca Ha and Julie Tseng, 2024-09-09
%   This file is part of MEGneto, see https://github.com/MabbottLab/MEGneto
%   for the documentation and details.
%
%% SETUP
    if exist('visitnum', 'var')
        this_output = [config.meta.project_path '/' config.meta.analysis_name '/' ...
            pid '/' sprintf('ses-%.2d', visitnum)]; % indicate subject-specific output folder path
    else
        this_output = [config.meta.project_path '/' config.meta.analysis_name '/' pid]; 
    end

    load([this_output '/out_struct.mat'])
    load([this_output '/step3_data_roi.mat'])
    
    if isstring(data_roi.label)
        data_roi.label = cellstr(data_roi.label);
    end

%% RUN DFC ANALYSIS   

    % set up parameters of dynamic FC analysis
    t_start   = config.step4c.toi(1); % start of time of interest
    t_end     = config.step4c.toi(2); % end of time of interest
    winsize   = config.step4c.winsize; 
    stepsize  = config.step4c.stepsize;
    
    %initializing output structure 
    connDFC           = [];
    connDFC.dimord    = 'chan_chan_freq';
    connDFC.label     = data_roi.label;
    connDFC.toi       = config.step4c.toi;
    connDFC.winsize   = winsize;
    connDFC.stepsize  = stepsize;
    connDFC.freqbands = config.step4c.freqbands;
    
    % send a warning if the step size is greater than the window size,
    % which would result in hugely overlapping windows
    if stepsize > winsize
        warning('Chosen step size is larger than window size');
    end 
    
    %%% FOR EACH FREQUENCY BAND --------------------------------------------
    for fq = 1:size(config.step4c.freqbands,1) 
        % filter data
        cfg           = [];
        cfg.bpfilter  = 'yes';
        cfg.bpfreq    = config.step4c.freqbands(fq,:);
        data_roi_filt = ft_preprocessing(cfg,data_roi);
        
        % defining variables for indexing and progress messages 
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
            cfg.latency       = [t t+winsize];
            data_roi_filt_toi = ft_selectdata(cfg,data_roi_filt); 
            
            % keep track of centre of windows 
            connDFC.win_centre{fq}(win_num) = t+winsize/2;
            % keep track of edges of windows 
            connDFC.win{fq}([1 2],win_num) = [t t+winsize];
            
            % Obtain power and cross-spectral density values necessary for
            % connectivity metric computation
            fprintf('Connectivity calculations for window: %.3f to %.3f s \n',t,t+winsize)
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

            % pull in labels which are impacted by reordering
            if ~isfield(conn, 'label') % if it hasn't already been done on a previous run
                connDFC.label = this_conn.label;
            end
            
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
    end   
    save([this_output '/step4c_connDFC.mat'], 'connDFC');     
end