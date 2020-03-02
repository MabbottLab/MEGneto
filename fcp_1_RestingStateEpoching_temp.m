function fcp_1_RestingStateEpoching_temp(paths)


%% SET UP LOGGING FILE

right_now = clock;
log_filename = [paths.conf_dir '/log_' sprintf('%d%d%d', right_now(1:3))];
diary(log_filename)

fprintf('\n\n%d:%d:%02.f       Now running **%s**.\n', ...
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
        %fcp1_output.trigger_figure  = 'triggerfigure.png';
    % data at various stages of cleaning
        fcp1_output.trial_cfg       = 'ft_meg_trl_cfg.json';
        fcp1_output.trial_cfgHM     = 'ft_meg_trl_cfgHM.json';
        fcp1_output.grad_cfg        = 'ft_meg_grad_cfg.json';
    % record keeping
        fcp1_output.subj_epochInfo  = 'subj_epoching_info.mat';
        fcp1_output.group_rmBadChan = 'group_rmBadChan.json';

%%% Preprocessing - Epoching -------------------------------------

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
    cfg             = ft_definetrial(cfg);
    



%%% HEAD MOTION CORRECTION ------------------------------------------------
    right_now = clock;
    fprintf('%d:%d:%02.f       Looking for excessive head motion...\n', right_now(4:6))    
	cfg.dataset = [paths.rawdata '/' subj_match.ds{ss}];
	cfg.trl = nx3 (matrix )
	cfg.channel                 = {'HLC0011','HLC0012','HLC0013', ...
        'HLC0021','HLC0022','HLC0023', ...
        'HLC0031','HLC0032','HLC0033'};
    headpos = ft_preprocessing(cfg); %sets a number of default pre-processing variables and configuration options using FT for the channels specified.
    fs = headpos.fsample;


    % find all HM values (in cm) in our segments and keep only the non-zero HM
    % values; this is with MEG turned off, but still recording
    HM = headpos.trial{1,1}(1,:);
    a = find(HM==0,1,'first');
    if isempty(a)
        new_HM=headpos.trial{1,1}(:,:); %process all data if there are no zero values
        time = 1:length(new_HM);
    else
        new_HM = headpos.trial{1,1}(:,1:a-1); %process all non-zero data if zero values exist
        time = 1:a-1;
    end


    try
        [~, ~, cfg, grad] = HeadMotionTool('Fieldtrip', cfg, ...
            'RejectThreshold', config.epoching.headMotion.thr, 'RejectTrials', true, 'CorrectInitial', true, ...
            'SavePictureFile', [paths.(subj_match.pid{ss}) '/' fcp1_output.fig_headmotion],...
            'GUI', false);
    catch
        warning('HeadMotionTool error!\n');
    end




%%% ARTIFACT DETECTION ----------------------------------------------------

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











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

