function Step4b_Connectivity(config, pid, visitnum)
%
% Step4b_Connectivity estimates functional connectivity (i.e., analyzes 
% the synchrony of signals from two regions) between regions of interest
% and over a given time period. 
%
% The most common functional connectivity metric used is "wpli_debiased"
% (see the ft_connectivityanalysis FieldTrip script for more details), but
% any of the FieldTrip ones are available for use. 
%
% INPUTS:
%   > config: 
%       struct, configured in your "main" script with
%       all analysis parameters and options and paths
%   > pid:
%       string, participant ID used to build the paths
%       to relevant data files and output folders
%   > visitnum:
%       int, visit number if longitudinal, optional arg)
%
% OUTPUTS:
%   > step4b_conn.mat:
%       FieldTrip style freq object with connectivity matrices in the
%       *.wpli_debiasedspctrm field. this field is a 4D matrix with the dim
%       order condition x roi x roi x frequency band
%
% Last updated by: Julie Tseng, 2024-09-09
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

load([this_output '/out_struct.mat']) % load metadata out struct
load([this_output '/step3_data_roi.mat']) % load data from step3_beamforming

% convert the label type from string to cellstring for fieldtrip function
% compatibility
if isstring(data_roi.label)
    data_roi.label = cellstr(data_roi.label);
end

%% RUN CONNECTIVITY ANALYSIS

%%% RUN CONNECTIVITY ANALYSIS ---------------------------------------------

    % instantiate "conn" variable to hold results
    conn = [];
    conn.dimord = 'chan_chan_freq';
    
    % handle multiple trial conditions and use first dimension of the
    % result to index the condition
    if ~isfield(config.step4b, 'conditions') % if no conditions
        config.step4b.conditions = {"all", 1:length(data_roi.trial)};
    else % if there are conditions, then find the trial indices associated to the condition
        config.step4b.conditions(:,2) = ...
            arrayfun(@(x) find(data_roi.trialinfo(:,1) == x & ... % find trial indices with condition x
            ~ismember(data_roi.trialinfo(:,2), ... % and 
            config.step4b.excludeTrialByIndex)), ... % exclude user-specified trials
            [config.step4b.conditions{:,2}], 'UniformOutput', false)';
    end

    % bookkeeping: keep track of the conditions in the final conn struct
    conn.conditions = config.step4b.conditions;
    
    for fq = 1:length(config.step4b.freqbands) % for each frequency band
        % filter signal to desired frequency band
        cfg             = [];
        cfg.bpfilter    = 'yes';
        cfg.bpfreq      = config.step4b.freqbands(fq,:);
        data_roi_filt   = ft_preprocessing(cfg, data_roi);

        for cond = 1:size(config.step4b.conditions, 1) % for each condition
              % use ft_select to select the timewindow of interest
              cfg           = [];
              cfg.latency   = config.step4b.toi; % timewindow of interest
              if isfield(config.step4b, 'conditions') % select trials corresponding to condition
                  cfg.trials = config.step4b.conditions{cond,2}; 
              end
              data_roi_filt_toi = ft_selectdata(cfg, data_roi_filt); % do the subselection

    %%% CALCULATE CONNECTIVITY ------------------------------------------------
                fprintf('Onto the connectivity calculations!\n')
                
                % first, grab the power/csd output which is necessary for
                % the connectivity analysis step
                cfg             = []; % set up config for connectivity calculation
                cfg.method      = 'mtmfft'; % fourier transform type
                cfg.output      = 'powandcsd'; % power and csd output
                cfg.channel     = 'all';
                cfg.trials      = 'all';
                cfg.keeptrials  = 'yes';
                cfg.taper       = 'hanning';
                cfg.foilim      = config.step4b.freqbands(fq,:);
                freq_filt       = ft_freqanalysis(cfg, data_roi_filt_toi); 

                % COMPUTE CONNECTIVITY METRIC
                cfg             = []; % set up config for computing connectivity metric
                cfg.method      = config.step4b.connmethod;
                this_conn       = ft_connectivityanalysis(cfg, freq_filt); % compute connectivity metric

                % RESHAPE INTO SOURCE X SOURCE CONN MAT AND STORE
                this_conn            = ft_checkdata(this_conn, 'cmbstyle', 'full');
                conn.freq(fq) = mean(this_conn.freq);
                
                % pull in labels which are impacted by reordering
                if ~isfield(conn, 'label') % if it hasn't already been done on a previous run
                    conn.label = this_conn.label;
                end

                % collapse frequencies within band
                if strcmp(config.step4b.bandavgmethod, 'max') % default to max across band
                    conn.(sprintf('%sspctrm',cfg.method))(cond,:,:,fq) = ...
                        squeeze(max(this_conn.(sprintf('%sspctrm', cfg.method)),[], 3));
                else % otherwise, use mean
                    conn.(sprintf('%sspctrm',cfg.method))(cond,:,:,fq) = ...
                        squeeze(mean(this_conn.(sprintf('%sspctrm', cfg.method)), 3));
                end
        end
    end

    save([this_output '/step4b_conn.mat'], 'conn');

end
