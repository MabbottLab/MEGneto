function fcp_1_TaskEpoching(paths)

% FCP_1_TASKEPOCHING will epoch MEG data into trials depending on the
% desired marker, detect and reject muscle/jump artifacts, and bad
% channels. This step will only reject trials with excessive head motion
% and muscle/jump artifacts. Bad channels are recorded but are repaired
% later on in the pipeline. 
% 
% NOTES:
%   - Ensure that subj_fcp1.csv is populated with the subject IDs of
%   included participants. 
%   - Desired parameters should be defined in the JSON config file.
%   User should double-check that the JSON config file is populated
%   appropriately, especially if a template JSON was copied over. 
%   - If the user wishes to browse the plotTriggers output, they must
%   set 'showFigure' to 'true; when the plotTriggers function is called.
%
% INPUTS:
%   paths               =   struct defining paths to data, participant
%                           folders, analysis folders, config files, etc. 
%
% OUTPUTS:
%   fcp1_output         = struct with locations of output files
%     .fig_headmotion   = 'headmotion.png': image of head motion graph
%     .trigger_figure   = 'triggerfigure.png': image of markers
%     .trial_cfg        = 'ft_meg_trl_cfg.json': fieldtrip cfg structure
%                          with the 'trl' parameter set to selected trials
%     .trial_cfgHM      = 'ft_meg_trl_cfgHM.json': head motion removed
%     .trial_grad_cfg   = 'ft_meg_grad_cfg.json': HM and 3rd order grad
%                          removed
%     .group_rmBadChan  = 'group_rmBadChan.json': lists of bad channs
%     .numtrls          = number of trials across participants
%     .HMremove_trls    = number of trials removed due to head motion
%     .noisy_trl        = muscle and jump artifact trial timestamps
%     .Nremove_trls     = total number of noisy trials removed
%     .bad_chann        = string array; bad channels detected for each 
%                         participant
%
% See also: DS_PID_MATCH, LOAD_PARTICIPANTS, WRITE_MATCH_IF_NOT_EMPTY, 
%           CHECK_CSV_HAS_EMPTY, PLOT_TRIGGERS, FT_READ_EVENT, 
%           FT_DEFINETRIAL, HEADMOTIONTOOL,FT_REJECTARTIFACT,
%           FT_ARTIFACT_MUSCLE, FT_ARTIFACT_JUMP, DETECTBADCHANNELS,
%           FT_PREPROCESSING

% Last updated by: Julie Tseng, 2020-01-07
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.

%% SET UP LOGGING FILE

right_now = clock;
log_filename = [paths.conf_dir '/log_' sprintf('%02.f%02.f%02.f', right_now(1:3))];
diary(log_filename)

fprintf('\n\n%02.f:%02.f:%02.f       Now running **%s**.\n', ...
    right_now(4:6), mfilename)

%% SETUP: LOAD CONFIG, PARTICIPANTS, CHECK FOR FULL DATASET, OUTPUTS

% load config JSON with analysis parameters
config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp1';

% check for matched MRI and MEG data
subj_match = ds_pid_match(paths,step);
if isempty(subj_match) % if there are no full sets of data
    error('No participants selected')
end

% save the final participants list to the subj_match_fcp1 CSV
write_match_if_not_empty(paths,step);

% check for multiple *.ds folders for participants
if length(unique(subj_match.pid)) ~= length(subj_match.pid)
    error('More than one ds per participant!')
end

% initialize output files
    % images
        fcp1_output.fig_headmotion  = 'headmotion.png';
        fcp1_output.trigger_figure  = 'triggerfigure.png';
    % data at various stages of cleaning
        fcp1_output.trial_cfg       = 'ft_meg_trl_cfg.json';
        fcp1_output.trial_cfgHM     = 'ft_meg_trl_cfgHM.json';
        fcp1_output.grad_cfg        = 'ft_meg_grad_cfg.json';
    % record keeping
        fcp1_output.subj_epochInfo  = 'subj_epoching_info.mat';
        fcp1_output.group_rmBadChan = 'group_rmBadChan.json';

%% EPOCHING
for ss = 1:length(subj_match.ds) % for each participant that has both MEG and MRI data
    
    fprintf('\n\n==================================\n...DS_FILE: %s\nSUBJECT: %s\n', ...
        subj_match.ds{ss}, subj_match.pid{ss});

%%% GRAB T0 MARKERS -------------------------------------------------------
    right_now = clock;
    fprintf('%02.f:%02.f:%02.f       Finding t0 markers...\n', right_now(4:6))

    numt0marker = plotTriggers(...
        [paths.rawdata '/' subj_match.ds{ss}], ...               % path to *.ds folder
        config.task.trialdef.markers.t0marker, 'savePath', ...   % consult config for t0 marker definition
        [paths.(subj_match.pid{ss}) '/' fcp1_output.trigger_figure], ... % save marker figure
        'showFigure', false); % set showFigure to true for debugging and viewing the plot 
    
    % if there were less than 5 markers found, throw a warning
    if numt0marker < 5
        warning('Not enough markers found - check marker file!')
        continue
    else % otherwise, display how many markers were found
        fprintf('\n %d markers were found for your specified t0markers\n ',...
            numt0marker)
    end
    
%%% EPOCHING --------------------------------------------------------------
    right_now = clock;
    fprintf('%02.f:%02.f:%02.f       Epoching into trials...\n', right_now(4:6))

    cfg             = []; % set up config parameters for ft_definetrial
    cfg.dataset     = [paths.rawdata '/' subj_match.ds{ss}]; 
    cfg.trialfun    = config.taskFunc; 
    cfg.trialdef    = config.task.trialdef;
    cfg.continuous  = 'yes';
    cfg             = ft_definetrial(cfg); 
    
    cfg_orig                    = cfg; % keep the original epoched data
    fcp1_output.numtrls{ss,1}   = length(cfg_orig.trl); % record num trials

%%% HEAD MOTION CORRECTION ------------------------------------------------
    right_now = clock;
    fprintf('%02.f:%02.f:%02.f       Looking for excessive head motion...\n', right_now(4:6))    

    try
        [~, ~, cfg, grad] = HeadMotionTool('Fieldtrip', cfg, ... % use the HeadMotionTool to display head movement info and remove bad trials
            'RejectThreshold', config.epoching.headMotion.thr, 'RejectTrials', true, 'CorrectInitial', true, ...
            'SavePictureFile', [paths.(subj_match.pid{ss}) '/' fcp1_output.fig_headmotion],...
            'GUI', false); 
    catch
        warning('HeadMotionTool error!\n');
    end
    
    % record number of trials removed due to head motion
    fcp1_output.HMremove_trls{ss,1} = length(cfg_orig.trl)-length(cfg.trl);
    % throw a warning if there are more than 90% of trials removed
    if fcp1_output.HMremove_trls{ss,1} > length(cfg_orig.trl)*0.9
        warning('\n\n \t\t Check head motion for this subject!\n\n')
        continue
    end
    
    save_to_json(cfg,...                % save the HM-removed data
        [paths.(subj_match.pid{ss}) '/' fcp1_output.trial_cfgHM],...
        true);

%%% ARTIFACT DETECTION ----------------------------------------------------
    right_now = clock;
    fprintf('%02.f:%02.f:%02.f       Detecting muscle and jump artifacts...\n', right_now(4:6))

    if config.cleaningOptions.artifact.detection == 1 % if user indicated they wish to perform artifact detection
        
        %%% Muscle Artifacts %%%
        % set up config for muscle artifact detection
        cfg.artfctdef.muscle.bpfilter    = config.cleaningOptions.artifact.muscle.bpfilter;
        cfg.artfctdef.muscle.bpfreq      = config.cleaningOptions.artifact.muscle.bpfreq;
        cfg.artfctdef.muscle.bpfiltord   = config.cleaningOptions.artifact.muscle.bpfiltord;
        cfg.artfctdef.muscle.bpfilttype  = config.cleaningOptions.artifact.muscle.bpfilttype;
        cfg.artfctdef.muscle.hilbert     = config.cleaningOptions.artifact.muscle.hilbert;
        cfg.artfctdef.muscle.boxcar      = config.cleaningOptions.artifact.muscle.boxcar;
        cfg.artfctdef.muscle.cutoff      = config.cleaningOptions.artifact.muscle.cutoff;
        cfg.artfctdef.muscle.trlpadding  = config.cleaningOptions.artifact.muscle.trlpadding; 
        cfg.artfctdef.muscle.fltpadding  = config.cleaningOptions.artifact.muscle.fltpadding;
        cfg.artfctdef.muscle.artpadding  = config.cleaningOptions.artifact.muscle.artpadding;
        [cfg, muscle_artifact]           = ft_artifact_muscle(cfg); % detect muscle artifacts
        
        %%%% Jump Artifacts %%%
        % set up config for jump artifact detection 
        cfg.artfctdef.jump.cutoff        = config.cleaningOptions.artifact.jump.cutoff;
        [cfg, jump_artifact]             = ft_artifact_jump(cfg); % detect jump artifacts
        
%%% ARTIFACT REJECTION ----------------------------------------------------
        % set up config for artifact rejection
        cfg.artfctdef.reject             = 'complete'; % remove complete trials
        cfg                              = ft_rejectartifact(cfg); % reject artifacts
        fcp1_output.noisy_trl{ss,1}      = [muscle_artifact; jump_artifact]; % record noisy trials
        fcp1_output.Nremove_trls{ss,1}   = length(cfg.trlold)-length(cfg.trl); % record number of removed trials
    end
    
    % save cleaned data progress
    save_to_json(cfg,... 
        [paths.(subj_match.pid{ss}) '/' fcp1_output.trial_cfg],...
        true);
    save_to_json(grad,...
        [paths.(subj_match.pid{ss}) '/' fcp1_output.grad_cfg],...
        true);

%%% BAD CHANNEL DETECTION -------------------------------------------------
    right_now = clock;
    fprintf('%02.f:%:02.f%02.f       Detecting bad channels...\n', right_now(4:6))
    
    cfg             = []; % set up config for detecting bad channels
    cfg.dataset     = [paths.rawdata '/' subj_match.ds{ss}];
    cfg.bchthr      = 60; % threshold; 75-85 quantile
    cfg.sections    = 3; % divide into 3 sections
    badChannels     = detectBadChannels(cfg,paths.name); % detect bad channels
    
    % record bad channels for that participant
    fcp1_output.bad_chann{ss,1} = badChannels;
    save_to_json(badChannels,...
        [paths.(subj_match.pid{ss}) '/badChannels.json'],...
        true);

%%% CELEBRATORY MESSAGE ---------------------------------------------------
    right_now = clock;
    fprintf('%02.f:%02.f:%02.f       Done subject %s!\n', ...
        right_now(4:6), subj_match.pid{ss})
    close all    
end % repeat for next participant

%%% FLAG PARTICIPANTS WITH MANY BAD CHANNELS ------------------------------

% who has more than 15 bad channels detected?
bad_subj = subj_match.pid(cellfun('length',fcp1_output.bad_chann) > 15); 
for i=1:length(bad_subj)
    warning([bad_subj{i},' has more than 15 BAD CHANNELS.']);
end

% save all output
save_to_json(fcp1_output, [paths.anout_grp '/fcp1_output.json'], true);

%% turn off diary
right_now = clock;
fprintf('%02.f:%02.f:%02.f       Done running **%s**.\n', ...
    right_now(4:6), mfilename)
diary off

%% let the users know that epoching is complete
sendEmail("epoching", string(config.contact));

end
