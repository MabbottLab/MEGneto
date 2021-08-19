%% Master run file
% Use the following code to run each step of the MEG pipeline in the
% presented sequence.

% Rename this file with something sensible like
% main_[analysis-name]_[date].m to denote your current analysis.

% Store this file somewhere sensible like in the project folder (whose 
% path is the "project_path", which is specified by the user in the code 
% block below), in the corresponding analysis folder (see "analysis_name"
% below), etc. 

%% Fill in the following paths for your project/analysis

% REPLACE WITH relevant folder paths
analysis_name             = '';	% folder where outputs will go
project_path              = ''; % folder where analysis folder will go (i.e., ~/project_folder/analysis/analysis_name/[config, analysis])
rawdata_path              = ''; % folder where raw MEG data is
mri_path                  = ''; % folder where MRI data is
megneto_path              = ''; % path to all MEGneto pipeline funtions
fieldtrip_path            = ''; % path to the FieldTrip toolbox/functions 
freeviewing_MEG_code      = ''; % path to MEG functions that had to be tailored to the freeviewing dataset
clipTimes_path            = ''; % path to folder where the output from functions used to extrapolate clip times
facesFeatureVector_path   = ''; % path to folder where the feature vector for faces v scenes, used in freeviewing_fcp_5_freqanalysis is used
OIRM_dataset_MRI          = ''; % path to folder where all OIRM MRI data is contained. This should be on the Carbon Drive under "OIRM_PILOT/MRI".

% add MEGneto and FieldTrip toolboxes to the path so that their contents
% can be accessed (without this, the functions cannot run)
addpath(freeviewing_MEG_code) % add path to freeviewing-edited functions
addpath(genpath(megneto_path)) % add path to remaining functions of the MEG repo that were not edited for freeviewing purposes 
addpath(fieldtrip_path) % add path to fieldtrip 
ft_defaults; % fieldtrip function to add correct sub-folders

% add paths to supplementary data
addpath(clipTimes); % add path to folder where the output from functions used to extrapolate clip times
addpath(facesFeatureVector_path); % add path to folder where the feature vector for faces v scenes, used in freeviewing_fcp_5_freqanalysis is used
addpath(OIRM_dataset_MRI); % add path to folder where all OIRM MRI data is contained. This should be on the Carbon Drive under "OIRM_PILOT/MRI".

%% if MATLAB crashes/closes and you need to re-load the paths variable, run the following:
% Note: the paths variable is the input to each fcp step

paths = loadjson(strcat(project_path, '/analysis/', analysis_name, '/config/paths.json'));

%%  fcp_0: setup (megne2setup)
%   This step involves defining paths to *.ds and *.mri data, as well as
%   output file/folder structure. 
%   **OUTPUTS**:
%       - paths struct:
%           Returned to MATLAB interface as a struct with address of each
%           folder (e.g., participant output, raw MRI data, etc.). This is
%           fed forward into each fcp_X function.
%       - [analysis_name].json: 
%           JSON file with all analysis parameters defined for every phase
%           of the pipeline. This can be copied from a template and then
%           modified by the user according to their analysis. 
%   **RESULTING FOLDER STRUCTURE**:
%       - project_folder (already exists)
%           - MRI_data_folder (already exists)
%           - MEG_data_folder (already exists)
%           - Analysis
%               - [analysis_name]
%                   - analysis (contains group & subject specific outputs)
%                       - group
%                       - subj_01
%                       - etc. (one folder per subject)
%                   - config (contains log files, csv files, JSON config file)

overwrite = false; % set to true if path set up for the same analysis should be re-done
paths = freeviewing_megne2setup(project_path, analysis_name, rawdata_path, mri_path, overwrite);

%% Interactive JSON config file
% Prior to running through the first step of the pipeline, users must
% define all relevant parameters in the JSON config file. To facilitate
% this process, the interactive_JSON_config function will guide users
% through the process of populating this guide. For more detail on the
% meaning of each parameter, please reference the "ConfigParams.md"
% document in the MEGneto repository on the Mabbott Lab GitHub.

interactive_JSON_config(paths, megneto_path) % run the interactive config function

%% subj_fcp1.csv population

% The following lines are used to auto-populate subj_fcp1.csv with
% participants who have .ds files. Should you wish to exclude certain
% participants for step 1 of the pipeline (fcp1), you can open the 
% subj_fcp1.csv file post auto-population and erase those participants.

MEG_ds = struct2table(dir(paths.rawdata)); % finds all MEG filenames
fid = fopen(paths.subj_fcp5, 'w'); % open subj_fcp1.csv
MEG_ds = MEG_ds.name(3:(height(MEG_ds))).'; % isolate only PIDs
fprintf(fid, '%s\n', MEG_ds{:}); % write each PID to file
fclose(fid); % close the file

%% fcp_1: task epoching, jump/muscle artifact detection, bad channel detection
% This step will epoch MEG data into trials based on desired markers,
% detect and reject trials with excessive head motion and muscle/jump
% artifacts, and detect and record (not reject) bad channels.

freeviewing_fcp_1_TaskEpoching(paths) % run first step of the MEG pipeline

% after fcp_1, check output for: not enough trials, excessive head motion, 
% and/or too many bad channels

%% fcp_2: ICA for artifact identification
% This step downsamples and filters/denoises MEG data, then carries out  
% ICA (if indicated in the config JSON file). 

% Remember to fill in participant IDs in the subj_fcp2.csv! Please navigate
% to subj_fcp1.csv, copy the list of participants, paste them in 
% subj_fcp2.csv and remove any participants you don't wish to include for 
% fcp2 and subsequent steps.

freeviewing_fcp_2_PreprocessingICA(paths)

%% fcp_2_5: human identifies bad ICA components, reject those components
% This interactive step guides the user to inspect ICA components from 
% fcp_2,and identify the ones that contain artifacts. After each 
% inspection, the user enters which components contain artifacts 

% Remember to fill in participant IDs in the subj_fcp2_5.csv! Please 
% navigate to subj_fcp2.csv, copy the list of participants, paste them in 
% subj_fcp2_5.csv and remove any participants you don't wish to include for 
% fcp2_5 and subsequent steps.

freeviewing_fcp_2_5_checkpoint(paths)

% remember to note participants to be excluded based on too many ICA artifact components!

%% fcp_3: bad channel repair
% This step repairs the bad channels that were detected in fcp_1 

% Remember to fill in participant IDs in the subj_fcp3.csv! Please navigate
% to subj_fcp2_5.csv, copy the list of participants, paste them in 
% subj_fcp3.csv and remove any participants you don't wish to include for 
% fcp3 and subsequent steps.

freeviewing_fcp_3_ChannelRepair(paths)

%% fcp_4: beamforming
% This step performs source reconstruction (the method is specified in the 
% JSON config file by the user) on the MEG data using anatomical
% MRI data. 

% Remember to fill in participant IDs in the subj_fcp4.csv! Please navigate
% to subj_fcp3.csv, copy the list of participants, paste them in 
% subj_fcp4.csv and remove any participants you don't wish to include for 
% fcp4 and subsequent steps.

% Uncomment the following to remove NAN channels for Run 1 data of
% participant P02. Previously this channel was found to be channel 182, but
% the user may change this channel number if they find otherwise.

% load([paths.OIRMP02 '/ft_meg_fullyProcessed.mat'],'-mat','data');
% cfg = [];
% cfg.channel = setdiff(1:length(data.label), 182);
% data = ft_selectdata(cfg, data);
% save([paths.OIRMP02 '/ft_meg_fullyProcessed'],'data','-v7.3')

freeviewing_fcp_4_beamforming(paths)

%% fcp_5: frequency analysis 
% This step performs frequency analysis to generate power spectrum data
% which can be used to test hypotheses based on spectral power.

% Remember to fill in participant IDs in the subj_fcp5.csv! Please navigate
% to subj_fcp4.csv, copy the list of participants, paste them in 
% subj_fcp5.csv and remove any participants you don't wish to include for 
% fcp5 and subsequent steps.

freeviewing_fcp_5_freqanalysis(paths);

%% The functions below are un-used in the freeviewing data analysis.
%% fcp_5: connectivity
% This step estimates functional connectivity (neural synchrony between
% regions) using a specified connectivity metric from the JSON config file.

fcp_5_taskconnectivity(paths);

%% Additional functions
% The sections below contain additional functions that may be run after
% the final fcp step of the pipeline. These functions serve as additional
% analysis/analysis preparation tools.

%% Preparation functions
% includes: make_NBS_ready, make_BNV_ready

%% make_NBS_ready 
% This function prepares a design matrix to serve as input to the
% NBS toolbox. The design matrix columns are the participant groups
% (e.g. "control", "surgery", etc.) and rows are participants. A "0" or
% "1" indicates whether the participant belongs to the group/column ("1") 
% or not ("0"). This function also outputs a list of participant
% (ordered_ppt_list.txt) in the analysis folder to inform the user of the
% design matrix's rows (same order as ppt list prestented in the .txt file).
% IMPORTANT: the by-frequency output files also have re-arranged
% participants that match the order of ordered_ppt_list.txt.

% Don't forget to include a nickname_ParticipantCategories.xlsx file 
% (where nickname corresponds to the subAnalysis_name variable below) in 
% your config folder. An example of this excel sheet with dummy variables
% is available in the templates folder (note that the column names of this
% file, in order, represent radiation, sugery, and typical development 
% controls). 
% Fill in the variables below which are the inputs to the function.

% Specify function inputs
group_names = NaN; % array of strings, e.g., ["surg", "rad", "control"], 
                    % exactly as it appears in the participant categories
                    % sheet
conn = NaN; % name of connectivity metric as a character array (must match 
%            the metric outlined in the file name of the connectivity
%            matrix .mat file). Can take on values including: 
%            "plv, "pli", "wpli", "wpli_debiased", "coh"

subAnalysis_name = NaN; % nickname for sub-analysis, e.g. 'TreatVsControl'
                        % this gets prepended to design_matrix and NBS data
                        % mat files. This should also match the name of
                        % the corresponding *_ParticipantCategories.xlsx

make_NBS_ready(paths, group_names, conn, subAnalysis_name)
%% make_BNV_ready
% This fuction creates *.node and *.edge files for viewing connectivity 
% results from PLS or NBS on BrainNet Viewer (BNV). 

% Don't forget to create the 'brainnet' struct containing nine user 
% specified parameters to pass into this function. 

% Specify function inputs
brainnet = NaN; % struct created by the user. Refer to documentation in the 
%                file 'make_BNV_ready.m' for a descrption of the parameters
%                in 'brainnet'.

make_BNV_ready(paths, brainnet)

%% Statistical analysis functions
% includes: bootTestDiffSeeds

%% bootTestDiffSeeds
% This function performs permutation-based significance testing (via 
% t-test or f-test using the max procedure) to build a null distribution 
% and control for Type 1 error.

% Note on T-max analysis extra input:
% If you wish to perform a T-max analsysis, you must include a .xlsx file
% with two groups that will be compared in the T-max analysis for 
% differences. You may use the make_NBS function to prepare this file
% where your input to make_NBS (the ParticipantCategories.xlsx) should have
% two groups you wish to include (e.g. "control" and "surg") and the
% participant IDs that fall under each group. 

% Specify function inputs
bootcfg.seed_regions = [1, 2, 3];      % numeric indices indicating the seed 
%                                       ROIs (e.g. if the AAL atlas is used, 
%                                       the default input [1, 2, 3] 
%                                       corresponds to the following regions 
%                                       ['left precentral gyrus', 'right 
%                                       precentral gyrus', 'left superior 
%                                       frontal gyrus, dorsolateral']. Note 
%                                       that for AAL atlas there are 90 
%                                       regions, so indices should take on 
%                                       values between 1-90). 
bootcfg.freq_band = 'gamma';          % frequency band of interst 
%                                       (e.g. 'alpha', 'beta', 'gamma', 
%                                       'theta')
bootcfg.num_bootstraps = 1000;        % number of desired bootstrap tests. 
%                                       Default is 1000.
bootcfg.thresh = 0.05;                % significance threshold for the 
%                                       p-value. Default is 0.05, can be 
%                                       altered to desired threshold 
%                                       by the user.
bootcfg.group_names = ["rad", "surg", "tdc"]; % array of strings, e.g., 
%                                       ["RAD", "SURG", "TDC"], 
%                                       exactly as they appear in the input
%                                       for group_names in the make_NBS
%                                       function.
bootcfg.nickname = 'nickname';         % subAnalysis name you set in make_NBS_ready

bootTestDiffSeeds(paths, bootcfg)

%% Summary functions
% includes: inspecting_results, getTrialSummary,
% getMarkerSummary,reactionTimes, trialSummary

%% inspecting_results
% This function allows the user to analyze the pipeline results by
% visualizing the data as specified by the user (in type.viz).

% Don't forget to include a ParticipantCategories.xlsx file in your
% paths.conf_dir folder and a 'type' struct as input.

% Specify function inputs
name = NaN; % name of the group output file the user wishes to inspect from 
%            fcp_5_freqanalysis (e.g. 'fcp_5_powspctrm_blcorrected.mat') or
%            fcp_5_taskconnectivity 
%            (e.g. 'fcp_5_allParticipants_conn_mats_wpli_debiased.mat')
type = NaN; % user-created struct with specifications for inspecting the 
%            data. Refer to the documentation in 'inspecting_results.m' 
%            for an example of the parameters included in the 'type'. 

inspecting_results(paths, name, type)

%% getTrialSummary
% This function summarizes trial information for each participant, such as
% the number of trials, the number of trials removed due to head motion, 
% number of trials removed due to noise, etc. 

% Specify function inputs
num_markers = NaN; % number of events expected 
                   % (total number of times stimulus is presented)
thresh = 25;       % percentage indicating what percentage of trials removed 
                   % is unacceptable. Here, 25 is the lab's convention.

getTrialSummary(paths, num_markers, thresh)

%% getMarkerSummary
% This function summarizes the markers for each participant and generates 
% a cluster which describes redundancies in the markers.

getMarkerSummary(paths)

%% reactionTimes
% This function extracts participant reaction times (time between stimulus 
% presentation and button press). Currently, this function is used for the 
% Speed of Thinking TP & EL dataset.

% This function is run by uncommenting certain lines within
% fcp_1_TaskEpoching.m. These lines include: 96, 97, 98, 102, 141, 142,
% 237, 238). Please read the reactionTimes.m file in the dev_functions
% folder for more details and information on where the output is stored. 

%% trialSummary
% This function generates a summary of the number of correct trials,
% incorrect trials and total trials per participant. This function must be
% run after at least fcp_1. Please read the reactionTimes.m file in the 
% dev_functions folder for more details and information on where the output 
% is stored. 

trialSummary(paths)
