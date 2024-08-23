%% MEGneto ANALYSIS TEMPLATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  The purpose of this "main" script is to act as a master pipeline script
%  and embedded location for you to note down quirks relevant to your
%  specific analysis. You should create a "main" script for each of your
%  analyses, as it will encode specific information about the markers
%  you're epoching around, threshold values, ROI indices, etc. 

%  To get started, copy this template into your
%  ../derivatives/meg_projectname_yourname/ folder, name it something clear
%  and obvious (e.g., main_motorRT.m), and start modifying below. 

%% Re-run if MATLAB has closed and you have already run setup before

cfg = []; % cfg is a struct that documents analysis options you've chosen and relevant filepaths

% in the "meta" field, encode relevant paths/metadata
    cfg.meta.analysis_name      = 'your_analysis';
    cfg.meta.project_path       = '/full/path/to/output/folder'; 
    cfg.meta.rawdata_path       = '/full/path/to/datasets/folder';
    cfg.meta.megneto_path       = '/full/path/to/MEGneto';
    cfg.meta.fieldtrip_path     = '/full/path/to/fieldtrip';
    cfg.meta.config_path        = [cfg.meta.project_path '/' ...
                                   cfg.meta.analysis_name ...
                                   '/config']; % auto-builds config folder path
    cfg.meta.participants_csv   = [cfg.meta.config_path '/participants_list.csv']; % path to list of participants
    
% add MEGneto and FieldTrip toolboxes to the path so that their contents
% can be accessed (without this, the functions cannot run)
    addpath(genpath(cfg.meta.megneto_path))
    addpath(cfg.meta.fieldtrip_path)
    ft_defaults; % fieldtrip function to add correct sub-folders

%% FIRST TIME SETUP ONLY: making directories and participants_list.csv

%  After you run this once to create directories, I suggest commenting out
%  this cell to avoid accidentally re-running it again

%  Check whether analysis folder exists and if not, create it
if ~exist([project_path '/' analysis_name], 'dir')
    mkdir([project_path '/' analysis_name])
    mkdir([project_path '/' analysis_name '/config'])
end

% Helper snippet of code to help people grab all the participant IDs available
% Manually edit as needed for your folder path / structure
ds_paths     = glob([cfg.meta.rawdata_path '/completion/pattern']);            % Use "glob" and some completion pattern (e.g., WMP/ses*/meg/*nback) to grab all participant IDs with MEG task data
out          = cellfun(@(x) strsplit(x, '/'), ds_paths, 'UniformOutput', false); % Split into participant ID and visit only
participants = table(string(cellfun(@(x) x{6}, out, 'UniformOutput', false)), ... % Convert into a table
                     strip(string(cellstr(ds_paths)),'right','/'), ...
                      string(cellfun(@(x) x{7}(end), out, 'UniformOutput', false)), ...
                     'VariableNames', {'ParticipantID', 'ds_path', 'VisitNum'}); % assumes VisitNum
writetable(participants, [cfg.meta.project_path '/' cfg.meta.analysis_name '/config/participants_list.csv']);
    
%% CONFIGURE YOUR ANALYSIS PARAMETERS HERE FOR EACH STEP %%%%%%%%%%%%%%%%%

%  Just as there is a cfg.meta field for the metadata, there is an
%  analogous field for each "Step" of the pipeline. Configure the options
%  for each step below. 

% Step1_Epoching %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % task definition
        cfg.step1.isTask = 1; % 1 = task, 0 = rest
        cfg.step1.maxRest = NaN; % NaN when isTask = 1
    % epoch definition
        % identify marker name to epoch around
        % use ft_read_event on the *.ds folder to probe list of markers if
        % you don't know which one you're going to use
        cfg.step1.trialdef.eventtype = 'MarkerName'; 
        
        % define how much pre/post marker to grab in seconds
        cfg.step1.trialdef.prestim = 1; % sec
        cfg.step1.trialdef.poststim = 1; % sec

        % which "trial function" to use for epoching?
        % you can develop your own custom one if you need to do more than
        % just epoch around some trial: for example, if you want to track
        % the trial index, or only include correct trials
        cfg.step1.trialdef.trialfun = 'my_trial_function'; 

        % head motion thresholding, see HeadMotionTool function
        cfg.step1.headMotionThr = 10; % head motion threshold in cm

        % do you want to check for muscle and jump artifacts?
        cfg.step1.clean_muscleJump = 1; % 1 = yes, 0 = no

% Step2_ICA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % do you want to run the ICA correction process? 
        cfg.step2.icaClean = 1; % 1 = yes, 0 = no
     
% Step3_Beamforming %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % dipole grid resolution, defines dipole spacing in source grid
        cfg.step3.templateRes = 6; % in mm

    % atlas details
        % is the atlas in the subject's native space (e.g., Glasser)?
        cfg.step3.nativeSpace = 1; % 1 = yes, 0 = no
        cfg.step3.atlas = 'mmp'; % provide atlas name

        cfg.step3.normLeadfield = 'no'; % set to 'yes' if analyzing resting state
        
    % source interpolation onto regions of interest
        % what ROI indices are you interested in? 
        cfg.step3.ROIs = [7, 163, 22, 5];

        % how do you want to combine signal from multiple dipoles belonging
        % to a particular region? PCA or mean?
        cfg.step3.combineDipoles = 'pca'; % or 'mean'
    
% Step4a_FrequencyAnalysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % specify the frequencies of interest (FOI)
        cfg.step4a.foi = start_fq:step:end_fq; % e.g., 2:2:100 = [2, 4, 6, ..., 100]
    
    % specify the frequency bands, particularly if you're only interested
    % in 1 or 2 of them
        cfg.step4a.freqbands = [1, 3; ...
                        4, 8; ...
                        8, 12; ...
                        13, 30; ...
                        30, 70; ...
                        70, 100];

    % IF YOU WANT TO RUN A POWER ANALYSIS (NO TIMEWINDOWS)
        cfg.step4a.type = 'power'; % power = overall power for the whole timewindow
        cfg.step4a.toi = [0.1 1]; % timewindow of interest in sec
        cfg.step4a.boi = [-1 -0.1]; % baseline window in sec
        cfg.step4a.baselinetype = 'absolute'; % options: absolute, relative, relchange
    % IF YOU WANT TO RUN A TIME WINDOW ANALYSIS (TFR)
        cfg.step4a.type = 'tfr'; % power at multiple timewindows
        cfg.step4a.toi = -1:0.01:1; % timewindows of interest in sec (e.g., -1 to 1 sec in increments of 0.01 sec) 
        cfg.step4a.boi = [-1 -0.1]; % baseline window in sec
        cfg.step4a.baselinetype = 'relchange'; % see ft_freqbaseline for options

% Step4b_Connectivity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cfg.step4b.connmethod = 'wpli_debiased';
    cfg.step4b.bandavgmethod = 'max';
    cfg.step4b.toi = [0 1]; 
    cfg.step4b.conditions = {"face_target", 1;
                            "thing_target", 2};
    cfg.step4b.excludeTrialByIndex = [1, 2, 3]; % indexes of trials to exclude
    cfg.step4b.freqbands = [1, 3; ...
                    4, 8; ...
                    8, 12; ...
                    13, 30; ...
                    30, 70; ...
                    70, 100];

% TO DO LIST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Group connectivity analysis 

    
%% actual commands

% run individually
Step1_Epoching(cfg, participants.ParticipantID{1}, participants.ds_path{1}, 1)
Step2_ICA(cfg, participants.ParticipantID{1}, "run", 1);

% run batch check on all participants
for p = participants.ParticipantsID
    Step2_ICA(cfg, p, "check", 1)
end

% regress bad components and fix bad channels
Step2_ICA(cfg, participants.ParticipantID{1}, "fix", 1);

% mark fiducial points on MRI
for p = participants.ParticipantsID
    Prep_T1(cfg, p)
end

% run beamforming
Step3_Beamforming(cfg, participants.ParticipantID{1}, 1)
    
% run power or tfr analysis
Step4a_FrequencyAnalysis(cfg, participants.ParticipantsID{1})

% run connectivity analysis
Step4b_Connectivity(cfg, participants.ParticipantID{1}, 1)
    
%% localizing oscillatory sources 

cfg                  = [];

