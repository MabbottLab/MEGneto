function fcp_1_RestingStateEpoching(paths)

% fcp_1_RestingStateEpoching will epoch MEG data into 2s trials, 
%  detect and reject muscle/jump artifacts, and 
% 
% NOTES:
%   - Ensure that subj_fcp1.csv is populated with the subject IDs of
%   included participants. 
%   - Desired parameters should be defined in the JSON config file.
%   User should double-check that the JSON config file is populated
%   appropriately, especially if a template JSON was copied over. 
%
% INPUTS:
%   paths               =   struct defining paths to data, participant
%                           folders, analysis folders, config files, etc. 
%
% OUTPUTS:
%   fcp1_output         = struct with locations of output files
%     .fig_headmotion   = 'headmotion.png': image of head motion graph
%     .trial_cfg        = 'ft_meg_trl_cfg.json': fieldtrip cfg structure
%                          with the 'trl' parameter set to selected trials
%     .trial_cfgHM      = 'ft_meg_trl_cfgHM.json': head motion removed
%     .trial_grad_cfg   = 'ft_meg_grad_cfg.json': HM and 3rd order grad
%                          applied
%     .group_rmBadChan  = 'group_rmBadChan.json': lists of bad channs
%     .numtrls          = number of trials across participants
%     .HMremove_trls    = number of trials removed due to head motion
%     .noisy_trl        = muscle and jump artifact trial timestamps
%     .Nremove_trls     = total number of noisy trials removed
%     .bad_chann        = string array; bad channels detected for each 
%                         participant
%
% See also: DS_PID_MATCH, WRITE_MATCH_IF_NOT_EMPTY, PLOTTRIGGERS, 
%           HEADMOTIONTOOL, DETECTBADCHANNELS

% Last updated by: Sonya Bells, 2020-02-13
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.

%% SET UP LOGGING FILE

right_now = clock;
log_filename = [paths.conf_dir '/log_' sprintf('%d%d%d', right_now(1:3))];
diary(log_filename)

fprintf('\n\n%d:%d:%02.f       Now running **%s**.\n', ...
    right_now(4:6), mfilename)

%% SETUP: LOAD CONFIG, PARTICIPANTS, CHECK FOR FULL DATASET, OUTPUTS

%% Epoching
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
        %fcp1_output.trigger_figure  = 'triggerfigure.png';
    % data at various stages of cleaning
        fcp1_output.trial_cfg       = 'ft_meg_trl_cfg.json';
        fcp1_output.trial_cfgHM     = 'ft_meg_trl_cfgHM.json';
        fcp1_output.grad_cfg        = 'ft_meg_grad_cfg.json';
    % record keeping
        fcp1_output.subj_epochInfo  = 'subj_epoching_info.mat';
        fcp1_output.group_rmBadChan = 'group_rmBadChan.json';

%% EPOCHING - looping through participants

for ss = 1:length(subj_match.ds) % for each participant
    
    fprintf('\n\n==================================\n...DS_FILE: %s\nSUBJECT: %s\n', ...
        subj_match.ds{ss}, subj_match.pid{ss});

	%%% EPOCHING Resting State-----------------------------------------------
    fprintf('Epoching into Segments...\n')

    % FT_REDEFINETRIAL allows you to adjust the time axis of your data, i.e. to change from stimulus-locked to response-locked.
    cfg 			= [];
    cfg.dataset     = [paths.rawdata '/' subj_match.ds{ss}]; 
    cfg.length		= 2; %epoch sections of 2s
    cfg.overlap 	= 0.5;
    cfg             = ft_redefinetrial(cfg,data);
    
    % This step removes  2) few segments of data at the end of the recording (which only contails 0s) 
    % 1) DC-component (later if fcp_2)
    cfg.trial 		= 1:(numel(cfg.trial)-6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cfg_orig                    = cfg; % keep the original epoched data
    fcp1_output.numtrls{ss,1}   = length(cfg_orig.trial); % record num trials

%%% HEAD MOTION CORRECTION ------------------------------------------------
    right_now = clock;
    fprintf('%d:%d:%02.f       Looking for excessive head motion...\n', right_now(4:6))    

    try
        [~, ~, cfg, grad] = HeadMotionTool('Fieldtrip', cfg, ...
            'RejectThreshold', config.epoching.headMotion.thr, 'RejectTrials', true, 'CorrectInitial', true, ...
            'SavePictureFile', [paths.(subj_match.pid{ss}) '/' fcp1_output.fig_headmotion],...
            'GUI', false);
    catch
        warning('HeadMotionTool error!\n');
    end
    
    % record number of trials removed due to head motion
    fcp1_output.HMremove_trls{ss,1} = length(cfg_orig.trial)-length(cfg.trial);
    % throw a warning if there are more than 90% of trials removed
    if fcp1_output.HMremove_trls{ss,1} > length(cfg_orig.trial)*0.9
        warning('\n\n \t\t Check head motion!!! \n\n')
        continue
    end
    
    save_to_json(cfg,...                % save the HM-removed data
        [paths.(subj_match.pid{ss}) '/' fcp1_output.trial_cfgHM],...
        true);


%%% ARTIFACT DETECTION ----------------------------------------------------
    right_now = clock;
    fprintf('%d:%d:%02.f       Detecting muscle and jump artifacts...\n', right_now(4:6))

    if config.cleaningOptions.artifact.detection == 1
        
        %%% Muscle Artifacts %%%
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
        [cfg, muscle_artifact]           = ft_artifact_muscle(cfg);
        
        %%%% Jump Artifacts %%%
        cfg.artfctdef.jump.cutoff        = config.cleaningOptions.artifact.jump.cutoff;
        [cfg, jump_artifact]             = ft_artifact_jump(cfg);
        
%%% ARTIFACT REJECTION ----------------------------------------------------
        cfg.artfctdef.reject             = 'complete'; % remove complete trials
        cfg                              = ft_rejectartifact(cfg);
        fcp1_output.noisy_trl{ss,1}      = [muscle_artifact; jump_artifact];
        fcp1_output.Nremove_trls{ss,1}   = length(cfg.trlold)-length(cfg.trial);
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
    fprintf('%d:%d:%02.f       Detecting bad channels...\n', right_now(4:6))
   
    cfg             = [];
    cfg.dataset     = [paths.rawdata '/' subj_match.ds{ss}];
    cfg.bchthr      = 60; % threshold; 75-85 quantile
    cfg.sections    = 3; % divide into 3 sections
    badChannels     = detectBadChannels(cfg,paths.name);
    
    % record bad channels for that participant
    fcp1_output.bad_chann{ss,1} = badChannels;
    save_to_json(badChannels,...
        [paths.(subj_match.pid{ss}) '/badChannels.json'],...
        true);

%%% CELEBRATORY MESSAGE ---------------------------------------------------
    fprintf('\nDone subject %s! \n',subj_match.pid{ss})
    close all    
end % repeat for next participant

%%% FLAG PARTICIPANTS WITH MANY BAD CHANNELS ------------------------------

% who has more than 15 bad channels detected?
bad_subj = subj_match.pid(cellfun('length',fcp1_output.bad_chann) > 15); 
for i=1:length(bad_subj)
    warning([bad_subj{i},' has more than 15 BAD CHANNELS.']);
end

% save all output
disp('Saving...');
save_to_json(fcp1_output, [paths.conf_dir '/fcp1_output.json'], true);
disp('Done FCP_1.');

%% turn off diary
right_now = clock;
fprintf('%d:%d:%02.f       Done running **%s**.\n', ...
    right_now(4:6), mfilename)
diary off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end



