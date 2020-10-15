function fcp_5_freqanalysis(paths)

% FCP_5_FREQANALYSIS ..... [UNDER CONSTRUCTION]  
%
% INPUTS:
%   paths               =   struct defining paths to data, participant
%                           folders, analysis folders, config files, etc. 
%
% OUTPUTS:

%
% See also: 
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

%% RUN CONNECTIVITY ANALYSIS
  

for ss = 1:length(subj_match.ds)
    right_now = clock;
    fprintf('%d:%d:%02.f       Working on subject %s!\n', ...
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
    data.time{tt} = time_info(1):(1/srate):time_info(2);
    %%% FOR EACH AAL NODE/SOURCE
        for kk = 1:num_sources
            data.trial{tt}(kk,:) = catmatrix(:,tt,kk);
        end
    end
    
    %%% Timewindow analysis
    fprintf('Onto the timewindow analysis!\n')
    cfg             = [];
    cfg.output      = 'pow';
    cfg.channel     = 'all';
    cfg.trials      = 'all';
    cfg.taper       = 'hanning';
    cfg.foi         = [2:2:100];        
    cfg.method = 'mtmconvol';
    cfg.t_ftimwin = 4 ./cfg.foi;
    cfg.toi         = -1.5:0.05:1.5;    
    freq = ft_freqanalysis(cfg, data);
    
    if ss == 1
        pow_spctrm = NaN(length(subj_match.ds), num_sources, length(cfg.foi), length(cfg.toi));
    end
    
    %%% baseline correct
    cfg                 = [];
    cfg.baseline        = [-1.5 -1];
    cfg.baselinetype    = 'relative';
    freq                = ft_freqbaseline(cfg, freq);

    pow_spctrm(ss,:,:,:) = freq.powspctrm;
end 

save([paths.anout_grp '/fcp_5_powspctrm_blcorrected.mat'],'pow_spctrm','-mat','-v7.3')


end