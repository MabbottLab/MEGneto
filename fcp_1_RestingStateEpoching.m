function fcp_1_RestingStateEpoching(paths)


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
% data at various stages of cleaning
fcp1_output.trial_cfg       = 'ft_meg_trl_cfg.json';
% fcp1_output.trial_cfgHM     = 'ft_meg_trl_cfgHM.json';
fcp1_output.grad_cfg        = 'ft_meg_grad_cfg.json';
% record keeping
fcp1_output.subj_epochInfo  = 'subj_epoching_info.mat';
fcp1_output.group_rmBadChan = 'group_rmBadChan.json';

%%% Preprocessing - Epoching -------------------------------------

for ss = 1:length(subj_match.ds) % for each participant
    
    fprintf('\n\n==================================\n...DS_FILE: %s\nSUBJECT: %s\n', ...
        subj_match.ds{ss}, subj_match.pid{ss});

    
%     cfg           = [];
%     cfg.dataset     = [paths.rawdata '/' subj_match.ds{ss}]; 
%     cfg.savemat     = 'no';
%     cfg.linefreq  = 60;
%     cfg.plotunit  = 300;
%     ft_qualitycheck(cfg)
%     
%%% EPOCHING Resting State-----------------------------------------------
    fprintf('Epoching into Segments...\n')
    
    % load raw data
    cfg                      = [];
    cfg.continuous           = 'yes';
    cfg.dataset              = [paths.rawdata '/' subj_match.ds{ss}]; 
    data                     = ft_preprocessing(cfg);
    
    % trim to remove recordings with no participant
    last_sample              = find(data.trial{1,1}(4,:) == 0, 1, 'first')-1;
    data.hdr.nSamples        = last_sample;
    data.time{1,1}           = data.time{1,1}(1:last_sample);
    data.trial{1,1}          = data.trial{1,1}(:,1:last_sample);
    data.sampleinfo(2)       = last_sample;
    
    % prep data struct for head motion detection
    cfg                      = [];
    cfg.continuous           = 'yes';
    cfg.dataset              = [paths.rawdata '/' subj_match.ds{ss}]; 
    cfg.trialfun             = config.taskFunc; % rest_trialfun.m
    cfg.trialdef.triallength = config.epoching.period;
    cfg.trialdef.overlap     = 0.5; % proportion overlap
    cfg.trialdef.endbound    = last_sample;
    data_epoched             = ft_definetrial(cfg); 
    
    % prep data struct for artifact detection
    cfg                      = [];
    cfg.length               = config.epoching.period;
    cfg.overlap              = 0.5;
    data_RAW_EPOCHED         = ft_redefinetrial(cfg, data);
     
    % record number of trials
    fcp1_output.numtrls{ss,1}   = length(data_epoched.trl); 

%%% ARTIFACT DETECTION ----------------------------------------------------

    if config.cleaningOptions.artifact.detection == 1
        
        %%% Muscle Artifacts %%%
        cfg                              = [];
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
        [cfg, muscle_artifact]           = ft_artifact_muscle(cfg, data_RAW_EPOCHED);
        
        %%%% Jump Artifacts %%%
        cfg.artfctdef.jump.cutoff        = config.cleaningOptions.artifact.jump.cutoff;
        [cfg, jump_artifact]             = ft_artifact_jump(cfg, data_RAW_EPOCHED);
        
%%% REJECTED ARTIFACT TRIALS ----------------------------------------------

        % use ft_rejectartifact to reject identified artifact trials
        cfg.artfctdef.reject             = 'complete'; % remove complete trials
        data_ARTFCT                      = ft_rejectartifact(cfg, data_RAW_EPOCHED);
        
        % take down which trials were rejected for later
        [~, rejected_artfct_trls]        = setdiff(data_ARTFCT.cfg.trlold(:,1), data_ARTFCT.cfg.trl(:,1));
        
        % record info on how many trials
        fcp1_output.noisy_trl{ss,1}      = [muscle_artifact; jump_artifact];
        fcp1_output.Nremove_trls{ss,1}   = length(rejected_artfct_trls);
    
 %%% HEAD MOTION CORRECTION ------------------------------------------------
        right_now = clock;
        fprintf('%d:%d:%02.f       Looking for excessive head motion...\n', right_now(4:6))      
        
        % use head motion tool on prepped data struct
        try        
            [~, ~, data_HMREMOVE, grad] = HeadMotionTool('Fieldtrip', data_epoched, ...
                'RejectThreshold', config.epoching.headMotion.thr, 'RejectTrials', true, 'CorrectInitial', true, ...
                'SavePictureFile', [paths.(subj_match.pid{ss}) '/' fcp1_output.fig_headmotion],...
                'GUI', false);
        catch
            warning('HeadMotionTool error!\n');
        end    
        
        % save gradiometer info
        save_to_json(grad, ...
            [paths.(subj_match.pid{ss}),'/' fcp1_output.grad_cfg], true);
        
        % take down which trials were rejected for later
        [~, rejected_HM_trls]        = setdiff(data_epoched.trl(:,1), data_HMREMOVE.trl(:,1));

        % record number of trials removed due to head motion       
        fcp1_output.HMremove_trls{ss,1} = fcp1_output.numtrls{ss,1}-length(data_HMREMOVE.trl);

        % throw a warning if there are more than 90% of trials removed
        if fcp1_output.HMremove_trls{ss,1} > fcp1_output.numtrls{ss,1}*0.9
            warning('\n\n \t\t Check head motion!!! \n\n')
            continue
        end
        
%%% REJECT ALL ARTIFACT/HEAD MOTION TRIALS AND DEMEAN ---------------------
        cfg                     = [];
        cfg.trials              = setdiff(1:fcp1_output.numtrls{ss,1}, unique([rejected_artfct_trls; rejected_HM_trls]))';
        cfg.demean              = 'yes';
        data                    = ft_preprocessing(cfg, data_RAW_EPOCHED);
        
        
%%% SAVE CLEANED, EPOCHED DATA (need to decide which format) --------------
        % save as *.mat
        save(fullfile(paths.(subj_match.pid{ss}),'/data_clean'), 'data');
        
        % save as *.json 
        save_to_json(data, ...
            [paths.(subj_match.pid{ss}),'/' fcp1_output.trial_cfg], true);
    end

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
end 

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

end


