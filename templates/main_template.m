%% Master run file
% Use the following to set up path names and run each step of the pipeline.
% Rename file with something sensible like main_[analysis-name]_[date].m
% Store in project folder

%% if you need to re-load paths in case MATLAB crashes/closes

% REPLACE WITH relevant folder paths
analysis_name = 'your_analysis_name';	% folder where outputs will go
project_path = 'your_project_folder'; % folder where analysis folder will go (i.e., ~/project_folder/analysis/analysis_name/[config, analysis])
rawdata_path = 'meg_path'; % folder where raw MEG data is
mri_path = 'mri_path'; % folder where MRI data is
megneto_path = 'megneto_path';
fieldtrip_path = 'fieldtrip_path

% add MEGneto and FieldTrip toolboxes
addpath(genpath(megneto_path))
addpath(fieldtrip_path)
ft_defaults % fieldtrip function to add correct sub-folders

% re-load your paths variable (input to each fcp step)
paths = loadjson(strcat(project_path, '/analysis/', analysis_name, '/config/paths.json');

%%  fcp_0: setup
%   This step involves defining paths to *.ds and *.mri data, as well as
%   output file/folder structure. 
%   **OUTPUTS**:
%       - paths struct:
%           Returned to MATLAB interface as a struct with address of each
%           folder (e.g., participant output, raw MRI data, etc.). This is
%           fed forward into each fcp_X function.
%       - [analysis_name].json: 
%           JSON file with all analysis parameters defined for every phase
%           of the pipeline. This can be copied from a template then
%           modified by the user according to their analysis. 
%   **RESULTING FOLDER STRUCTURE**:
%       - project_folder (already exists)
%           - MRI_data_folder (already exists)
%           - MEG_data_folder (already exists)
%           - Analysis
%               - [analysis_name]
%                   - analysis
%                       - group
%                       - subj_01
%                       - 
%                   - config
overwrite = false;
paths = megne2setup(project_path, analysis_name, rawdata_path, mri_path, overwrite);

%% fcp_1: task epoching, jump/muscle artifact detection, bad channel detection

% uncomment the following two lines if you'd like to auto-populate subj_fcp1.csv
% MEG_ds = struct2table(dir(paths.rawdata)); % finds all MEG filenames
% writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp1')); % writes them to file

fcp_1_TaskEpoching(paths) % run first step

% remember to check output for: not enough trials, excessive head motion, too many bad channels

%% fcp_2: ICA for artifact identification

% remember to fill in participant IDs in the subj_fcp2.csv!
fcp_2_PreprocessingICA(paths)

%% fcp_2_5: human identifies bad ICA components, reject those components

fcp_2_5_checkpoint(paths)

% remember to note participants to be excluded based on too many ICA artifact components!

%% fcp_3: bad channel repair

% remember to fill in participant IDs in the subj_fcp3.csv!
fcp_3_ChannelRepair(paths)

%% fcp_4: beamforming
% remember to fill in participant IDs in the subj_fcp4.csv!
fcp_4_beamforming(paths)

%% fcp_5: frequency analysis 
% remember to fill in participant IDs in the subj_fcp5.csv!
fcp_5_freqanalysis(paths);

%% fcp_5: connectivity
fcp_5_taskconnectivity(paths);


