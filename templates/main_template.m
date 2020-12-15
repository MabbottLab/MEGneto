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
analysis_name = 'your_analysis_name';	% folder where outputs will go
project_path = 'your_project_folder'; % folder where analysis folder will go (i.e., ~/project_folder/analysis/analysis_name/[config, analysis])
rawdata_path = 'meg_path'; % folder where raw MEG data is
mri_path = 'mri_path'; % folder where MRI data is
megneto_path = 'megneto_path'; % path to all MEGneto pipeline funtions
fieldtrip_path = 'fieldtrip_path'; % path to the FieldTrip toolbox/functions 

% add MEGneto and FieldTrip toolboxes to the path so that their contents
% can be accessed (without this, the functions cannot run)
addpath(genpath(megneto_path))
addpath(fieldtrip_path)
ft_defaults; % fieldtrip function to add correct sub-folders

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
paths = megne2setup(project_path, analysis_name, rawdata_path, mri_path, overwrite);

%% fcp_1: task epoching, jump/muscle artifact detection, bad channel detection
% This step will epoch MEG data into trials based on desired markers,
% detect and reject trials with excessive head motion and muscle/jump
% artifacts, and detect and record (not reject) bad channels.

% uncomment the following lines if you'd like to auto-populate subj_fcp1.csv
% MEG_ds = struct2table(dir(paths.rawdata)); % finds all MEG filenames
% fid = fopen(paths.subj_fcp1, 'w'); % open subj_fcp1.csv
% MEG_ds = MEG_ds.name(3:(height(MEG_ds))).'; % isolate only PIDs
% fprintf(fid, '%s\n', MEG_ds{:}); % write each PID to file
% fclose(fid) % close the file

fcp_1_TaskEpoching(paths) % run first step of the MEG pipeline

% after fcp_1, check output for: not enough trials, excessive head motion, too many bad channels

%% fcp_2: ICA for artifact identification
% This step downsamples and filters/denoises MEG data, then carries out  
% ICA (if indicated in the config JSON file). 

% remember to fill in participant IDs in the subj_fcp2.csv! If you would
% like to auto-populate subj_fcp2.csv (i.e. include all participants), copy
% the commented lines of code from fcp_1 and chaged "path.subj_fcp1" to
% "paths.subj_fcp2"

fcp_2_PreprocessingICA(paths)

%% fcp_2_5: human identifies bad ICA components, reject those components
% This interactive step guides the user to inspect ICA components from 
% fcp_2,and identify the ones that contain artifacts. After each 
% inspection, the user enters which components contain artifacts 

% remember to fill in participant IDs in the subj_fcp2.csv! If you would
% like to auto-populate subj_fcp2.csv (i.e. include all participants), copy
% the commented lines of code from fcp_1 and chaged "path.subj_fcp1" to
% "paths.subj_fcp2_5"

fcp_2_5_checkpoint(paths)

% remember to note participants to be excluded based on too many ICA artifact components!

%% fcp_3: bad channel repair
% This step repairs the bad channels that were detected in fcp_1 

% remember to fill in participant IDs in the subj_fcp3.csv!If you would
% like to auto-populate subj_fcp2.csv (i.e. include all participants), copy
% the commented lines of code from fcp_1 and chaged "path.subj_fcp1" to
% "paths.subj_fcp3"

fcp_3_ChannelRepair(paths)

%% fcp_4: beamforming
% This step performs source reconstruction (the method is specified in the 
% JSON config file by the user) on the MEG data using anatomical
% MRI data. 

% remember to fill in participant IDs in the subj_fcp4.csv!If you would
% like to auto-populate subj_fcp2.csv (i.e. include all participants), copy
% the commented lines of code from fcp_1 and chaged "path.subj_fcp1" to
% "paths.subj_fcp4"

fcp_4_beamforming(paths)

%% fcp_5: frequency analysis 
% This step performs frequency analysis to generate power spectrum data
% which can be used to test hypotheses based on spectral power.

% remember to fill in participant IDs in the subj_fcp5.csv!If you would
% like to auto-populate subj_fcp2.csv (i.e. include all participants), copy
% the commented lines of code from fcp_1 and chaged "path.subj_fcp1" to
% "paths.subj_fcp5"

fcp_5_freqanalysis(paths);

%% fcp_5: connectivity
% This step estimates functional connectivity (neural synchrony between
% regions) using a specified connectivity metric from the JSON config file.

fcp_5_taskconnectivity(paths);

%% Additional functions
% The section(s) below contain additional functions that may be run after
% the final fcp step of the pipeline. These functions serve as additional
% analysis/analysis preparation tools.

%% make_NBS_ready 
% This function prepares a design matrix to serve as input to the Matlab
% NBS toolbox. The design matrix columns are the participant groups
% (e.g. "control", "surgery", etc.) and rows are participants. A "0" or
% "1" indicates whether the participant belongs to the group/column ("1") 
% or not("0"). 

% Don't forget to include a ParticipantCategories.xlsx file in your
% paths.conf_dir folder. Fill in the variables below which are input to the
% function.

% Specify inputs to the function
group_names = []; % array of strings, e.g., ["surg", "rad", "control"], 
%                   exactly as they appear in folder names 
conn = ""; % name of connectivity metric as a character array (must match 
%            the metric outlined in the file name of the connectivity
%            matrix .mat file, e.g. "wpli_debiased")

make_NBS_ready(paths, group_names, conn)