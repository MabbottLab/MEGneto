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

% Note that this step uses freeviewingTaskTrialFunc where the output of the
% clip times extraction step (done prior to using this freeviewing MEGneto
% pipeline), is loaded to facilitate trial epoching.

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

% Note that this step saves two parts of beamforming data for each run of
% participant data. This is done to facilitate beamforming at a finer grid
% resolution. 

% Uncomment the following to remove NAN channels for Run 1 data of
% participant P02. Previously this channel was found to be channel 182, but
% the user may change this channel number if they find otherwise.

% load([paths.OIRMP02 '/ft_meg_fullyProcessed.mat'],'-mat','data');
% cfg = [];
% cfg.channel = setdiff(1:length(data.label), 182);
% data = ft_selectdata(cfg, data);
% save([paths.OIRMP02 '/ft_meg_fullyProcessed'],'data','-v7.3')

freeviewing_fcp_4_beamforming(paths)

%% fcp_4_5: reinserting trials
% This step stacks the output of the beamforming step for each participant 
% (i.e., it combines all beamforming from both runs of a participant's MEG 
% data) and inserts blank columns for trials that were rejected in the task 
% epoching step. 

% This step is not meant to be run as a function. Instead, the user should
% open the script and run each section of code as specified at the top of
% that file.

%% fcp_5: frequency analysis 
% This step performs frequency analysis to generate power spectrum data
% which can be used to test hypotheses based on spectral power.

% Remember to fill in participant IDs in the subj_fcp5.csv! Please navigate
% to subj_fcp4.csv, copy the list of participants, paste them in 
% subj_fcp5.csv and remove any participants you don't wish to include for 
% fcp5 and subsequent steps.

% Note that this step loads in a feature vector that discriminates
% categories of visual content. The user must specify which category
% should be analyzed, as only one category is run at a time. For instance,
% if the user wishes to analyse the faces category first, they must specify
% that in the freeviewing_fcp_5_freqanalysis function and alter the name of
% the saved file, and then run the function. Then, they must change the 
% category to scenes, change the output name to reflect the scenes
% category, and re-run the function. If the user wishes to analyse all
% categories at once (i.e., not separate it out into groups), the user must
% specify that as well. All of this is specified in lines 135-137.

freeviewing_fcp_5_freqanalysis(paths);

%% fcp_5_5: statistical analysis
% This step performs statistical analysis and creates visual 
% representations of results for each aim and hypothesis of the 
% OIRM Faces v Scenes Summer 2021 research project. 

% This step is meant to be run in sections, not as a function. The user 
% should open the script and refer to it for detailed instructions on
% each section of the script.