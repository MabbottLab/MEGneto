function fcp_2_PreprocessingICA(paths)

% FCP_2_PREPROCESSINGICA will prepare epoched data for ICA preprocessing,
% then carry out the ICA omitting signal from bad channels. This step also
% accounts for 3rd order gradients. 
% 
% NOTES:
%   - Check participants who had excessive head motion or excessive numbers
%   of bad channels. 
%   - Ensure that subj_fcp2.csv is populated with the subject IDs of
%   participants you want to include after checking over initial results. 
%   - Need to remove bad channels from ica - if not you will get complex 
%   numbers. Because during repair channels procedure bad channels are 
%   repaired according to neighbours, thus the new ones are not unique (no 
%   independent components).
%
% INPUTS:
%   paths               =   struct defining paths to data, participant
%                           folders, analysis folders, config files, etc. 
%
% OUTPUTS:
%   fcp2_output         = struct with locations of output files
%     .data_noisecorr   = mat file with gradients accounted for
%     .preprocessedData_cfg = ft configuration
%     .ICAcomp_cfg      = JSON of ICA components
%     .data_icacomp     = mat file of ICA components

%
% See also: 

% Last updated by: Julie Tseng, 2020-01-08
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.

%% SET UP LOGGING FILE

right_now = clock;
log_filename = [paths.conf_dir '/log_' sprintf('%d%d%d', right_now(1:3))];
diary(log_filename)

fprintf('\n\n%d:%d:%02.f       Now running **%s**.\n', ...
    right_now(4:6), mfilename)

%% SETUP: LOAD CONFIG, CHECK PPTS, FCP_1 OUTPUT

% load config JSON with analysis parameters
config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp2';

% check for matched MRI and MEG data
subj_match = ds_pid_match(paths,step);
ssSubjPath = @(x) paths.(subj_match.pid{x});
if isempty(subj_match) % if there are no full sets of data
    error('No participants selected')
end

% save the final participants list to the subj_match_fcp1 CSV
write_match_if_not_empty(paths,step);

% check for multiple *.ds folders for participants
if length(unique(subj_match.pid)) ~= length(subj_match.pid)
    error('More than one ds per participant!')
end

% load output from previous step
fcp1_output     = loadjson([paths.anout_grp '/fcp1_output.json']);

% initialize output files
    % data
        fcp2_output.preprocessedData_cfg  =  'ft_meg_data_cfg.mat';
        fcp2_output.data_noisecorr        =  'data_noisecorr.mat';
    % ica components
        fcp2_output.ICAcomp_cfg           =  'ft_icacomp.json'; % save ica components
        fcp2_output.data_icacomp          =  'icacomponents.mat';

%% PROCESSING

% determine which subjects to process
rangeOFsubj = 1:height(subj_match);
for ss = rangeOFsubj %1:rangeOFsubj

    fprintf('\n\n==================================\nSUBJECT: %s\n', subj_match.pid{ss});

    % check whether preprocessed files already exist
    if ~exist([ssSubjPath(ss) '/' fcp2_output.data_noisecorr], 'file')
        right_now = clock;
        fprintf('%d:%d:%02.f       Starting preprocessing...\n', right_now(4:6))

    %%% LOAD DATA -------------------------------------------------------------
        % load trial definition
        if exist([ssSubjPath(ss) '/' fcp1_output.trial_cfg], 'file')
            % if RESTING STATE
            if config.task.isRest == false 
                samples = loadjson([ssSubjPath(ss) '/' fcp1_output.trial_cfg]);
                cfg     = [];
                cfg.trl = samples.trl;
                disp('Epoched file found!');
                disp([ssSubjPath(ss) '/' fcp1_output.trial_cfg]);
            % if TASK
            else 
                % data where only excessive head motion removed
                if config.cleaningOptions.artifact.rmNoisyTrls == 0
                    samples = loadjson([ssSubjPath '/' fcp1_output.trial_cfgHM]);
                    cfg     = [];
                    cfg.trl = samples.trl;
                    disp('Epoched file found!');
                    disp([ssSubjPath(ss) '/' fcp1_output.trial_cfgHM]);
                % data where artifacts were also rejected
                elseif config.cleaningOptions.artifact.rmNoisyTrls == 1
                    samples = loadjson([ssSubjPath(ss) '/' fcp1_output.trial_cfg]);
                    cfg     = [];
                    cfg.trl = samples.trl;
                    disp('Epoched file found!');
                    disp([ssSubjPath(ss) '/' fcp1_output.trial_cfg]);
                end
            end
        else
            error('Epoched mat %s was not found', [ssSubjPath(ss) '/' fcp1_output.trial_cfg]);
        end

        % load data
        cfg.dataset = [paths.rawdata '/' subj_match.ds{ss}];

    %%% FILTER DATA -----------------------------------------------------------
        right_now = clock;
        fprintf('%d:%d:%02.f       Filtering data...\n', right_now(4:6))

        cfg.channel     = config.filteringParameters.channel;
        cfg.dftfilter   = config.filteringParameters.dftfilter;
        cfg.dftfreq     = config.filteringParameters.dftfreq;
        cfg.bpfilter    = config.filteringParameters.bpfilter;
        cfg.bpfreq      = config.filteringParameters.bpfreq;
        cfg.bpfiltord   = config.filteringParameters.bpfiltord;
        cfg.continuous  = 'yes';
        dataFiltered = ft_preprocessing(cfg);

    %%% ACCOUNTING FOR GRADIOMETERS -------------------------------------------
        if exist([ssSubjPath(ss) '/' fcp1_output.grad_cfg] , 'file')
            right_now = clock;
            fprintf('%d:%d:%02.f       Loading gradiometer config from file...\n', right_now(4:6))

            dataFiltered.grad       = loadjson([ssSubjPath(ss) '/' fcp1_output.grad_cfg]);
            dataFiltered.hdr.grad   = dataFiltered.grad;
        else
            warning('No modified gradiometer definitions found!');
        end

        % synthetic 3rd order grads - for noise reduction in CTF no ref channels stored
        cfg             = [];
        cfg.gradient    = 'G3BR';
        data_noisecorr  = ft_denoise_synthetic(cfg,dataFiltered);
        
        % Resample all datasets
        cfg             = [];
        cfg.resamplefs  = config.filteringParameters.sampleRate;
        cfg.detrend     = 'no';
        data_resamp     = ft_resampledata(cfg, data_noisecorr);
        data_noisecorr  = data_resamp;
        clear data_resamp;

        % save data with noise reduction thru 3rd order gradients
        save([ssSubjPath(ss) '/' fcp2_output.data_noisecorr], 'data_noisecorr', '-v7.3')
    else
        fprintf('Preprocessed data already exists! Skipping...\n');
    end
end

% remove participants
pid_fcp1    = readtable(paths.subj_fcp1_match); % get fcp_1 ppts
pid_fcp1    = pid_fcp1.Var1; % in a struct format for comparison
pid_fcp2    = subj_match.pid; % get fcp_2 ppts
match_func  = cellfun(@(x) ismember(x, pid_fcp2), pid_fcp1, 'UniformOutput', 0);
included    = find(cell2mat(match_func)); % get indices of included

for ss = rangeOFsubj
    
    % load filtered data
    load([ssSubjPath(ss) '/' fcp2_output.data_noisecorr]);

    % if you don't want to run ICA
    if config.cleaningOptions.artifact.icaClean == 0
        right_now = clock;
        fprintf('%d:%d:%02.f       No ICA requested, saving data...\n', right_now(4:6))

        % save data file
        data = data_noisecorr;
        clear data_noisecorr;
        save([ssSubjPath(ss) '/' fcp2_output.preprocessedData_cfg],'data','-v7.3');
        disp('Done.');
%%% ICA -------------------------------------------------------------------
    elseif config.cleaningOptions.artifact.icaClean == 1
        right_now = clock;
        fprintf('%d:%d:%02.f       Running ICA...\n', right_now(4:6))

        % downsample the data to speed up the next step
        cfg             = [];
        cfg.resamplefs  = 300;
        cfg.detrend     = 'no';
        data_resampICA  = ft_resampledata(cfg, data_noisecorr);

        % first, handle bad channels if they exist
        % load info about bad channels
        channel_check = dir([paths.(subj_match.pid{ss}) '/badChannels.json']);

        if channel_check.bytes > 5 % 5 = NULL, anything larger means populated
            badChanns   = loadjson([paths.(subj_match.pid{ss}) '/badChannels.json']); % get list of bad channels from fcp_1
            cfg         = [];
            cfg.channel = setdiff(data_noisecorr.label, badChanns);         % specify included channels
            data_resampICArmBadCh = ft_selectdata(cfg, data_resampICA);     % select channel data

            % Run ICA
            cfg          = [];
            cfg.channel  = 'MEG';
            cfg.method   = 'fastica'; % default and uses the implementation from EEGLAB
            comp         = ft_componentanalysis(cfg, data_resampICArmBadCh);
        else % if there are no bad channels, proceed
            cfg          = [];
            cfg.channel  = 'MEG';
            cfg.method   = 'fastica'; 
            comp         = ft_componentanalysis(cfg, data_resampICA);
        end
        % save the ICA components
        save([ssSubjPath(ss) '/' fcp2_output.data_icacomp],'comp','-v7.3')
    end

    close all
    
    right_now = clock;
    fprintf('%d:%d:%02.f       Done subject %s!\n', right_now(4:6), subj_match.pid{ss})

end

%%% SAVE fcp_2 output -----------------------------------------------------
save_to_json(fcp2_output,[paths.anout_grp '/fcp2_output.json'])

%% turn off diary
right_now = clock;
fprintf('%d:%d:%02.f       Done running **%s**.\n', ...
    right_now(4:6), mfilename)
diary off

%% let the users know
sendEmail("ICA", string(config.contact));


end
