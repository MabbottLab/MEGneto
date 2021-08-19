function freeviewing_fcp_5_freqanalysis(paths)

% Compared to the standard MEG pipeline, the following lines were edited:
% 83-85
% 117-128
% 135-137
% 147
% 158
% 161-170

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
% freq_structure        a structure for each ppt containing the power 
%                       spectrum and power spectrum's separate components
%
% OTHER: 
%    The user is required to load in a feature vector file that indicates
%    how visual stimuli was separate (see lines 81-83).

%
% See also: FT_FREQANALYSIS
%
% Last updated by: Julie Tseng, 2020-07-02
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.
%
%% SET UP LOGGING FILE

right_now = clock;
log_filename = [paths.conf_dir '/log_' sprintf('%02.f%02.f%02.f', right_now(1:3))];
diary(log_filename)

fprintf('\n\n%02.f:%02.f:%02.f       Now running **%s**.\n', ...
    right_now(4:6), mfilename)

%% SETUP

%%% PARTICIPANT IDS -------------------------------------------------------
% load config JSON with analysis parameters
config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp5';

% check is user indicated that time frequency analysis should be performed
if ~config.freqanalysis.include
    error("Your input to the JSON config file indicates you do not wish to run fcp_5_freqanalysis. If you wish to run this function, please change the config.freqanalysis.include field to 1.")
end

% check for matched MRI and MEG data
subj_match  = ds_pid_match(paths,step);
ssSubjPath  = @(x) paths.(subj_match.pid{x});

%% RUN FREQUENCY ANALYSIS

for ss = 1:length(subj_match.ds)
    right_now = clock;
    fprintf('%02.f:%02.f:%02.f       Working on subject %s!\n', ...
        right_now(4:6), subj_match.pid{ss})

%%% LOAD VIRTUAL SENSOR DATA ---------------------------------------------
    try
        load([ssSubjPath(ss) '/revamped_beamforming_results.mat'], '-mat'); 
    catch
        load([ssSubjPath(ss) '/AAL_beamforming_results.mat'], '-mat');
    end
    
    % load feature vector to separate clips by category of visual content
    facesFeatureVector = load([facesFeatureVector_path 'facesFeatureVector.mat']);
    facesFeatureVector = facesFeatureVector.';
    
    % define some dimensions
    num_samples = size(catmatrix, 1);
    num_trials  = size(catmatrix, 2);
    
    % collapsed across ROIs if indicated
    ROIs = config.freqanalysis.ROIs;
    if ~isempty(ROIs{1})
        catmatrix_collapsed = cellfun(@(x) nanmean(catmatrix(:,:,x), 3), ROIs, 'UniformOutput', false);
        catmatrix = cat(3, catmatrix_collapsed{:}); clear catmatrix_collapsed;
    end
    num_sources = size(catmatrix, 3);

    data    = [];
    time_info = config.task.trialdef.parameters.tEpoch;
    for src = 1:num_sources
        data.label{src} = sprintf('ROI%d', src);
    end

    %%% FOR EACH TRIAL
    for tt = 1:num_trials
    fprintf('Processing trial %d...\n', tt);
    data.time{tt} = time_info(1):(1/srate):time_info(2); % split the time (necessary for subsequent steps)
    %%% FOR EACH OF THE NODES/SOURCES
        for kk = 1:num_sources
            data.trial{tt}(kk,:) = catmatrix(:,tt,kk); % reformat source space data
        end
    end
    
    index = [];
    
    %%% Identify all NaN trials
    for trial = 1:length(data.trial)
        disp(sum(sum(isnan(data.trial{trial}))))
        if all(all(isnan(data.trial{trial})))
            index = [index trial];
        end
    end
    
    %%% Remove trials that are fully NaN
    data.trial(index) = [];
    data.time(index) = [];
    facesFeatureVector(index) = [];
    
    %%% Timewindow analysis
    fprintf('Onto the timewindow analysis!\n')
    cfg             = []; % set up config for the frequency analysis
    cfg.output      = 'pow';
    cfg.channel     = 'all';
    %%% uncomment line 122 and comment out line 123 if you are separating trials based on visual stimuli
    %cfg.trials      = find(facesFeatureVector == 1); % 1 for faces, 0 for other
    cfg.trials      = 'all';
    cfg.taper       = 'hanning';
    cfg.foi         = config.freqanalysis.foi; % frequencies of interest
    cfg.method      = config.freqanalysis.method;
    cfg.t_ftimwin   = 4 ./cfg.foi;
    cfg.toi         = config.freqanalysis.toi; % time of interest
    freq            = ft_freqanalysis(cfg, data); % perform frequency analysis
    
    if ss == 1
        pow_spctrm = NaN(length(subj_match.ds), num_sources, length(cfg.foi), length(cfg.toi)); % set up blank matrix to store power spectrum data
        freq_structure = cell(length(subj_match.pid), 1);
    end
    
    %%% baseline correct
    % necessary to control for general/random spikes in power
    cfg                 = []; % set up config for baseline correction
    cfg.baseline        = config.freqanalysis.baseline;
    cfg.baselinetype    = config.freqanalysis.baseline_type;
    freq                = ft_freqbaseline(cfg, freq); % perform baseline correction

    pow_spctrm(ss,:,:,:) = freq.powspctrm; % store the participant's power spectrum data
    freq_structure{ss,1} = freq; % store the participants structure
end 

%%% save the output: output file names are changed based on the input to
%%% the function. For instance, if line 121 is used to find faces, outputs
%%% should be named to reflect that. Similarly, if the TOI in the JSON
%%% config file is altered (e.g., to occur in 33msec increments) the file
%%% name should also reflect that.

%save([paths.anout_grp '/fcp_5_powspctrm_blcorrected_faces.mat'],'pow_spctrm','-mat','-v7.3') 
%save([paths.anout_grp '/fcp_5_freq_structure_faces.mat'], 'freq_structure', '-mat', '-v7.3')
save([paths.anout_grp '/fcp_5_powspctrm_blcorrected_33msecBin.mat'],'pow_spctrm','-mat','-v7.3') 
save([paths.anout_grp '/fcp_5_freq_structure_33msecBin.mat'],'freq_structure','-mat','-v7.3') 

%% turn off diary
right_now = clock;
fprintf('%d:%d:%02.f ============== Finished Processing ====================\n', ...
    right_now(4:6))
diary off

%% let the users know that connectivity analysis is complete
sendEmail("connectivity", string(config.contact));

end