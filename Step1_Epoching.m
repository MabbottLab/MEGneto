function Step1_Epoching(config, pid, ds_path)
%
% some sort of documentation will be here
% 1. epoch
% 2. artifact: head motion, muscle, jump
% ... excess cleanup if rest
% 3. bad channel detection
%
%% SETUP

this_output = [config.meta.project_path '/' config.meta.analysis_name '/' pid];
if ~exist(this_output, 'dir') 
    mkdir(this_output)
end

out = [];

%% EPOCHING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if config.step1.isTask % if this is a task
    %%% PLOT T0 MARKERS -------------------------------------------------------
        numt0marker = plotTriggers(...
            ds_path, ...               % path to *.ds folder without last '/'
            config.step1.trialdef.eventtype, ...
            'savePath', [this_output '/plot_markers.png'], ... % save marker figure
            'showFigure', false); % set showFigure to true for debugging and viewing the plot 
        
        % save the total number of markers
        out.step1.numTrials_epoched = numt0marker;

    %%% TRIAL SETUP --------------------------------------------------------------
        cfg             = []; % set up config parameters for ft_definetrial
        cfg.datafile    = ds_path;
        cfg.trialdef    = config.step1.trialdef;
        cfg             = ft_definetrial(cfg); 
    else % if this is a resting state scan     
        cfg                      = [];
        cfg.dataset              = ds_path;
        cfg.continuous           = 'yes';
        data                     = ft_preprocessing(cfg);

        % trim to remove recordings with no participant
        last_sample              = find(data.trial{1,1}(4,:) == 0, 1, 'first')-1;
        if ~isempty(last_sample)
            data.hdr.nSamples        = last_sample;
            data.time{1,1}           = data.time{1,1}(1:last_sample);
            data.trial{1,1}          = data.trial{1,1}(:,1:last_sample);
            data.sampleinfo(2)       = last_sample;
        else
            last_sample = data.sampleinfo(2);
        end
        
        % prep data struct for head motion detection
        cfg                      = [];
        cfg.continuous           = 'yes';
        cfg.dataset              = [paths.rawdata '/' subj_match.ds{ss}]; 
        cfg.origFs               = data.hdr.Fs;
        cfg.trialfun             = config.taskFunc; % rest_trialfun.m
        cfg.trialdef.triallength = config.epoching.period;
        cfg.trialdef.overlap     = config.epoching.overlap; % proportion overlap
        cfg.trialdef.endbound    = last_sample;
        cfg                      = ft_definetrial(cfg); 
        clear data last_sample
    end % end of epoching

%% ARTIFACT HANDLING: HEAD MOTION, MUSCLE, JUMP %%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% head motion correction --------------------------------------------
    try
        [~, ~, cfg, grad] = HeadMotionTool('Fieldtrip', cfg, ... % use the HeadMotionTool to display head movement info and remove bad trials
            'RejectThreshold', config.step1.headMotionThr, ... % threshold in mm
            'RejectTrials', true, ... % remove trials with excessive motion
            'CorrectInitial', true, ... % use different definition for initial head loc
            'SavePictureFile', [this_output '/headmotion.png'], ... % save figure
            'GUI', false); 
    catch
        warning('HeadMotionTool error!\n');
    end
    
    % track the number of trials with rejected head motion
    out.step1.numTrials_headMotionRejected = ...
        out.step1.numTrials_epoched - length(cfg.trl);

    %%% muscle and jump artifacts -----------------------------------------
    if config.step1.clean_muscleJump == 1 % if user indicated they wish to perform artifact detection
        
        %%% Muscle Artifacts %%%
        % use default muscle detection configuration
        cfg.artfctdef.muscle.bpfilter    = 'yes';
        cfg.artfctdef.muscle.bpfreq      = [110, 140]; % frequency range
        cfg.artfctdef.muscle.bpfiltord   = 8; % filter order
        cfg.artfctdef.muscle.bpfilttype  = 'but'; % butterworth filter
        cfg.artfctdef.muscle.hilbert     = 'yes';
        cfg.artfctdef.muscle.boxcar      = 0.2;
        cfg.artfctdef.muscle.cutoff      = 30;
        if config.step1.isTask % if task
            cfg.artfctdef.muscle.trlpadding  = 0.5;
            cfg.artfctdef.muscle.fltpadding  = 0.1;
            cfg.artfctdef.muscle.artpadding  = 0.1;
        else % if rest
            cfg.artfctdef.muscle.trlpadding  = 0; 
            cfg.artfctdef.muscle.fltpadding  = 0;
            cfg.artfctdef.muscle.artpadding  = 0;
        end
        [cfg, muscle_artifact]           = ft_artifact_muscle(cfg); % detect muscle artifacts
        
        %%%% Jump Artifacts %%%
        % set up config for jump artifact detection 
        cfg.artfctdef.jump.cutoff        = 35;
        [cfg, jump_artifact]             = ft_artifact_jump(cfg); % detect jump artifacts
        
        % save these artifact settings
        out.step1.artfctdef              = cfg.artfctdef;
        
    %%% actual rejection of trials ----------------------------------------
        % set up config for artifact rejection
        cfg.artfctdef.reject             = 'complete'; % remove complete trials
        cfg                              = ft_rejectartifact(cfg); % reject artifacts
        
        % do some bookkeeping
        out.step1.artfctdef.artsamples   = [muscle_artifact; jump_artifact];
        out.step1.numTrials_artRejected  = numt0marker - length(cfg.trl);
        out.step1.numTrials_final        = length(cfg.trl);
    end

%%% EXCESS TRIAL CLEAN-UP (REST)-------------------------------------------
    % if you have resting state data and specified a max length of time
    if ~config.step1.isTask && ~isnan(config.step1.maxRest) 
        ntrials                          = size(cfg.trl, 1); % grab number of trials
        maxSamp                          = config.step1.maxRest*cfg.origFs; % grab max sample
        lastTrial                        = find(cfg.trl(:,2) <= maxSamp, 1, 'last'); % find trial containing end sample
        
        if ntrials > lastTrial % trim if we have more than we need
            cfg.trl = cfg.trl(1:lastTrial,:);
        end
    end

    % save this trial definition
    out.step1.cfg = cfg;
    out.step1.grad = grad;

%%% BAD CHANNEL DETECTION -------------------------------------------------
    cfg             = []; % set up config for detecting bad channels
    cfg.dataset     = ds_path;
    cfg.bchthr      = 60; % threshold; 75-85 quantile
    cfg.sections    = 3; % divide into 3 sections
    if config.step1.isTask
        badChannels     = detectBadChannels(cfg,'Task'); % detect bad channels
    else
        badChannels     = detectBadChannels(cfg,'Rest');
    end
    
    % save the bad channel detection parameters
    out.step1.badChanDef.bcthr = cfg.bchthr;
    out.step1.badChanDef.sect  = 3;
    out.step1.badChanDef.out   = badChannels;

%% PREPROCESS THE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cfg             = out.step1.cfg; % call up again our trial def
    cfg.channel     = {'MEG', 'MEGREF', 'REFGRAD', 'REFMAG'};
    % line noise removal
    cfg.dftfilter   = 'yes'; % line noise removal using discrete fourier transform? 
    cfg.dftfreq     = [60, 120]; % line noise frequencies
    % overall signal filtering
    cfg.bpfilter    = 'yes'; % bandpass filter = specify range
    cfg.bpfreq      = [1, 150]; % bandpass range delta to high gamma
    cfg.bpfiltord   = 5; % filter parameter
    data_filtered   = ft_preprocessing(cfg);
    
    % incorporate gradiometers signal
    data_filtered.hdr.grad = grad;
    cfg = [];
    cfg.gradient = 'G3BR';
    data_noisecorr = ft_denoise_synthetic(cfg, data_filtered);
    
    % resample to lower hz (freq res = 1/2 sampling, so 300Hz default)
    cfg = [];
    cfg.resamplefs = 300;
    cfg.detrend = 'no';
    data_clean = ft_resampledata(cfg, data_noisecorr);
    
    % save some necessary outputs
    out.step1.filt.linefreq = [60, 120];
    out.step1.filt.bandpass = [1, 150];
    out.step1.filt.bandpassOrder = 5;
    out.step1.resamplefs = cfg.resamplefs;
    
    % save out struct
    save([this_output '/out_struct.mat'], 'out');
    
    % save current state of data
    save([this_output '/step1_data_clean.mat'], 'data_clean');
    
end
