function [ ] = fcp_1_TaskEpoching(paths)
% fcp_1taskepoching
% Version 3 Ming Scott
% 2019 June 19
% Version 2 Sonya Bells
% 2016 November 1
% Version 1 Simeon Wong, Anne Keller
% 2016 March 7
%
% INPUT FILES:
% - raw .ds (MEG 'RAW' dataset with at least 3 min of continuous recording)
%
% OUTPUT FILES:
% - ft_cfg.mat containing a fieldtrip cfg structure with the 'trl' parameter set to the selected
%   trials.
%
% Before running this function, make sure you have subj_fcp1.csv in the
% configs set up, and your preferred parameters in the analysis json

%% Read in config settings and match *.ds files to subject IDs
config = load_config(paths, paths.name);
config = config.config;
step = 'fcp1';
subj_ds = load_participants(paths,step);
pids = readtable(paths.all_subj_pids);
[subj_match, failure] = ds_pid_match(paths,step);

if isempty(subj_match)
    error('No participants selected')
end

write_match_if_not_empty(paths,step);

if length(unique(subj_match.pid)) ~= length(subj_match.pid)
    error('More than one ds per participant!')
end

%%% OUTPUTS %%% 
%%% Note: Only set up on inital run, edit config in config files afterwards
fcp1_output.fig_headmotion = 'headmotion.png';
fcp1_output.trial_cfg = 'ft_meg_trl_cfg.json';
fcp1_output.trial_cfgHM = 'ft_meg_trl_cfgHM.json';
fcp1_output.grad_cfg = 'ft_meg_grad_cfg.json';
fcp1_output.subj_epochInfo  = 'subj_epoching_info.mat';

fcp1_output.group_rmBadChan = 'group_rmBadChan.json';

% Trial timeseries figures
fcp1_output.trigger_figure = 'triggerfigure.png';
fcp1_output.timeseries_figure = @(trial) ['Trial_', sprintf('%04d',trial), '.png'];

save_to_json(fcp1_output, [paths.conf_dir '/fcp1_output.json'])

%% Run epoching
% determine which subjects to process
rangeOFsubj = 1:length(subj_match.ds);

for ss = rangeOFsubj
    
    %     if ~p.subj.include(ss) % skip subject if excluded from analysis
    %         continue
    %     end
    
    fprintf('\n\n==================================\nDS_FILE: %s\nSUBJECT: %s\n', ...
        subj_match.ds{ss}, subj_match.pid{ss});
    
    
    %%% Plot triggers %%%
%     numt0marker = plotTriggers(...
%         [paths.rawdata '/' subj_match.ds{ss}],...
%         config.task.trialdef.markers.t0marker, 'savePath',...
%         [paths.(subj_match.pid{ss}) '/' fcp1_output.trigger_figure],...
%         'showFigure', true); %for Debugging
    numt0marker = plotTriggers(...
        [paths.rawdata '/' subj_match.ds{ss}],... % for some reason, mine was a full file name already
        config.task.trialdef.markers.t0marker, 'savePath',...
        [paths.(subj_match.pid{ss}) '/' fcp1_output.trigger_figure],...
        'showFigure', true); % set showFigure to true for debugging
    if numt0marker < 5
        warning('Not enough markers found - check marker file!')
        continue
    else
        fprintf('\n %d markers were found for %s \n ',...
            numt0marker,config.task.trialdef.markers.t0marker)
    end
    
    %%% Check if there are missing Triggers %%%
    %Define structure of error checking
    %     p.trialdef.Error.includeError = 0;
    %     p.trialdef.Error.includeOnceError = 0;
    %     p.trialdef.Error.t0marker = 0;
    %     p.trialdef.Error2.t0marker = 0;
    %
    %     p = checkForMissingTriggers( p.subj.subj_ds{ss} , p);
    %     if ~p.trialdef.Error.t0marker(ss)
    %           continue
    %     else
    %% Do Epoching %%%
    
    cfg = [];
    cfg.dataset = [paths.rawdata '/' subj_match.ds{ss}];
    cfg.trialfun = config.taskFunc; %@taskTrialFun;
    cfg.trialdef = config.task.trialdef;
    cfg.continuous         = 'yes';
    cfg  = ft_definetrial(cfg);
    cfg_orig = cfg;
    %%% Head motion correction %%%
    % call headmotion tool
%     try
        [~, ~, cfg, grad] = HeadMotionTool('Fieldtrip', cfg, ...
            'RejectThreshold', config.epoching.headMotion.thr, 'RejectTrials', true, 'CorrectInitial', true, ...
            'SavePictureFile', [paths.(subj_match.pid{ss}) '\' fcp1_output.fig_headmotion],...
            'GUI', false);
%     catch
%         warning('HeadMotionTool error!');
%     end
    fcp1_output.Numtrls{ss,1} = length(cfg_orig.trl);
    fcp1_output.HMremove_trls{ss,1} = length(cfg_orig.trl)-length(cfg.trl);
    if fcp1_output.HMremove_trls{ss,1} > length(cfg_orig.trl)*0.9
        warning(sprintf('\n\n \t\t Check head motion!!! \n\n'))
        continue
    end
    
    
    cfg_HM = cfg;
    save_to_json(cfg,...
        [paths.(subj_match.pid{ss}) '/' fcp1_output.trial_cfgHM],...
        true);
    %    save(p.paths.trial_cfgHM((ss)), '-struct', 'cfg');
    %%
    if config.cleaningOptions.artifact.detection == 1
        %%%% find muscle artifacts %%%
        %Optimal for identifying muscle artifacts
        cfg.artfctdef.muscle.bpfilter    = 'yes';
        cfg.artfctdef.muscle.bpfreq      = [110 140];
        cfg.artfctdef.muscle.bpfiltord   = 8;
        cfg.artfctdef.muscle.bpfilttype  = 'but';
        cfg.artfctdef.muscle.hilbert     = 'yes';
        cfg.artfctdef.muscle.boxcar      = 0.2;
        
        cfg.artfctdef.muscle.cutoff      = 30; % default is 8 (makes enough epochs for good DS)
        cfg.artfctdef.muscle.trlpadding  = 0.5; % get errors in ft_artifact_muscle without these
        cfg.artfctdef.muscle.fltpadding  = 0.1;
        cfg.artfctdef.muscle.artpadding  = 0.1;
        [cfg,muscle_artifact] = ft_artifact_muscle(cfg);
        
        %%%% find jump artifacts %%%
        cfg.artfctdef.jump.cutoff        = 35 ; % default is 22 (makes enough epochs for good DS)
        [cfg,jump_artifact] = ft_artifact_jump(cfg);
        
        cfg.artfctdef.reject = 'complete'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
        %         cfg.artfctdef.eog.artifact = artifact_EOG; %
        %         cfg.artfctdef.jump.artifact = artifact_jump;
        %         cfg.artfctdef.muscle.artifact = artifact_muscle;
        cfg = ft_rejectartifact(cfg);
        total_artifacts = [muscle_artifact; jump_artifact];
        % add number of noise trials
        fcp1_output.noisy_trl{ss,1} = total_artifacts;
        fcp1_output.Nremove_trls{ss,1} = length(cfg.trlold)-length(cfg.trl);
    end
    
    %%% Save fieldtrip configuration %%%
    save_to_json(cfg,...
        [paths.(subj_match.pid{ss}) '/' fcp1_output.trial_cfg],...
        true);
    save_to_json(grad,...
        [paths.(subj_match.pid{ss}) '/' fcp1_output.grad_cfg],...
        true);
%     save(p.paths.trial_cfg((ss)), '-struct', 'cfg');
%     save(p.paths.grad_cfg((ss)), '-struct', 'grad');
    
    %%% Finding Bad Channels %%%
    disp('   ')
    disp('Detecting bad channels ...')
    disp('   ')
    cfg                         = [];
    cfg.dataset                 = [paths.rawdata '/' subj_match.ds{ss}];
    cfg.bchthr = 60;% 75-85 quantile
    cfg.sections = 3; % divids into 3 sections
    [badChannels,resMat] = detectBadChannels(cfg,paths.name);
    
    % add bad channels to p-structure
    fcp1_output.bad_chann{ss,1} = badChannels;
    fprintf('\nDone subject %s! \n',subj_match.pid{ss})
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all
    
end 
disp('   ');
bad_subj = subj_match.pid(cellfun('length',fcp1_output.bad_chann) > 15) ; % find subjects with more than 15 bad channels (alert user)
for i=1:length(bad_subj)
    warning([bad_subj{i},' has more than 15 BAD CHANNELS.']);
end


% save all bad channels (just in case)
all_bad_chann = fcp1_output.bad_chann;
save_to_json(all_bad_chann, [paths.anout_grp '/' fcp1_output.group_rmBadChan], true);

% save all output
save_to_json(fcp1_output, [paths.anout_grp '/fcp1_output'], true);
% save p-structure
% disp('Saving...');
% save(p.paths.p_strct,'p','-mat','-v7');
% disp('Done.');
end
