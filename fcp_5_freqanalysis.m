function fcp_5_freqanalysis(paths)

% FCP_5_FREQANALYSIS uses spectral analysis on time-frequency 
% representations of data to test hypotheses based on spectral power. The 
% virtual sensor data from the beamforming step is loaded in and frequency
% analysis is performed on various timewindows of the data (time window 
% analysis). Thus, for each subject and each interpolated atlas region, a 
% power spectrum is calculated and corrected to a baseline to control for 
% general/random spikes in power.
% This function is useful for cases where hypothesis are based on power in
% some general region of interest.
%
% INPUTS:
%   paths               =   struct defining paths to data, participant
%                           folders, analysis folders, config files, etc. 
%
% OUTPUTS:
% pow_spctrm:           a 4-D matrix containing power spectrum data, where
%                       dimensions are: 
%                       [participants] x [regions] x [frequency] x [time].

%
% See also: FT_FREQANALYSIS
%
% Last updated by: Julie Tseng, 2020-07-02
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.
%
%% SETUP

%%% PARTICIPANT IDS -------------------------------------------------------
% load config JSON with analysis parameters
config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp5';

% check for matched MRI and MEG data
subj_match  = ds_pid_match(paths,step);
ssSubjPath  = @(x) paths.(subj_match.pid{x});

%% RUN FREQUENCY ANALYSIS

for ss = 1:length(subj_match.ds)
    right_now = clock;
    fprintf('%02.f:%02.f:%02.f       Working on subject %s!\n', ...
        right_now(4:6), subj_match.pid{ss})

%%% LOAD VIRTUAL SENSOR DATA ----------------------------------------------
    load([ssSubjPath(ss) '/AAL_beamforming_results'], '-mat'); 
    
    % define some dimensions
    num_samples = size(catmatrix, 1);
    num_trials  = size(catmatrix, 2);
    num_sources = size(catmatrix, 3);
    
    data    = [];
    time_info = config.task.trialdef.parameters.tEpoch;
    for src = 1:num_sources
        data.label{src} = sprintf('AAL%d', src);
    end

    %%% FOR EACH TRIAL
    for tt = 1:num_trials
    fprintf('Processing trial %d...\n', tt);
    data.time{tt} = time_info(1):(1/srate):time_info(2); % split the time (necessary for subsequent steps)
    %%% FOR EACH AAL NODE/SOURCE
        for kk = 1:num_sources
            data.trial{tt}(kk,:) = catmatrix(:,tt,kk); % reformat source space data
        end
    end
    
    %%% Timewindow analysis
    fprintf('Onto the timewindow analysis!\n')
    cfg             = []; % set up config for the frequency analysis
    cfg.output      = 'pow';
    cfg.channel     = 'all';
    cfg.trials      = 'all';
    cfg.taper       = 'hanning';
    cfg.foi         = [2:2:100];        
    cfg.method = 'mtmconvol';
    cfg.t_ftimwin = 4 ./cfg.foi;
    cfg.toi         = -1.5:0.05:1.5; % time of interest    
    freq = ft_freqanalysis(cfg, data); % perform frequency analysis
    
    if ss == 1
        pow_spctrm = NaN(length(subj_match.ds), num_sources, length(cfg.foi), length(cfg.toi)); % set up blank matrix to store power spectrum data
    end
    
    %%% baseline correct
    % necessary to control for general/random spikes in power
    cfg                 = []; % set up config for baseline correction
    cfg.baseline        = [-1.5 -1];
    cfg.baselinetype    = 'relative';
    freq                = ft_freqbaseline(cfg, freq); % perform baseline correction

    pow_spctrm(ss,:,:,:) = freq.powspctrm; % store the participant's power spectrum data
end 

%%% save the output
save([paths.anout_grp '/fcp_5_powspctrm_blcorrected.mat'],'pow_spctrm','-mat','-v7.3') 


end