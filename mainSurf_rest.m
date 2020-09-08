%% Master run file - RESTING STATE PIPELINE

% Use the following to set up path names and run each step of the pipeline.
addpath(genpath('/home/megneto/MEGneto'))
addpath('/home/megneto/fieldtrip')
ft_defaults

% replace filepath with the one specific to your study analysis
paths = loadjson('/home/megneto/JT_TEST/analysis/rest/config/paths.json');

% FreeSurfer surfaces filepath
surf_path = '/home/megneto/data/sbells/meg_freesurfer/processed_freesurfer';

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

analysis_name       = 'rest';
project_path        = '/home/megneto/JT_TEST';
rawdata_path        = '/home/megneto/data/sbells/meg_freesurfer/MEG';
mri_path            = '/home/megneto/data/sbells/meg_freesurfer/MRIs';
paths               = megne2setup(project_path, analysis_name, rawdata_path, mri_path, false);

%% fcp_1: task epoching, jump/muscle artifact detection, bad channel det

% write participant IDs to the appropriate *.csv, helper lines
MEG_ds = struct2table(dir(paths.rawdata)); % get all possible IDs
writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp1')); % write them

% run analysis
fcp_1_RestingStateEpoching(paths)

%% fcp_2: ICA

% write PIDs to appropriate *.csv, ASSUMING no change in participants to be
% analyzed
writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp2'));

% run analysis
fcp_2_RestPreprocessingICA(paths)

%% fcp_2_5: INTERACTIVE - identify bad ICA components (human)

% write PIDs then modify according to who you want to remove
writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp2_5'));

% run interactive analysis
fcp_2_5_checkpoint(paths)

%% fcp_3: bad channel repair

% write PIDs then modify according to who you want to remove
writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp3'));

% run channel repair
fcp_3_ChannelRepair(paths)

%% fcp_4a: INTERACTIVE - register freesurfer surfaces to ctf mri
%  user interaction needed for ac-pc labeling and freesurfer determine
%  coordinate system (RAS)

% write PIDs for analysis
writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp4'));

% run registration
fcp_4a_register_CTF_Surf(paths,surf_path)

%%% SONYA'S NOTES ON OIRM PARTICIPANTS-------------------------------------
% no acpc:
%   good : 02,03,07,24,26
%   bad  : 06,11,15,16,17,27,28,

% w acpc (spm12)
%   good : 02,03,06,07,11,15,16,17,24,26,27,28,29,30,31,39,41,43,44,47
%   16??17??

%% fcp_4: beamforming surfaces

% write PIDs
writecell(MEG_ds.name(3:4), paths.('subj_fcp4'));

% run beamforming
fcp_4_Beamforming_Freq_Surf(paths)

% NOTE TO USERS: 
% freq_band : [min, max], center, +/- value
% -------------------------------------------
% theta     : [ 4,  8],   5,  2
% alpha     : [ 8, 13],  10,  2 
% beta      : [13, 30],  20,  7
% gamma     : [30, 48],  39,  9

%% fcp_5: connectivity

% write PIDs
writecell(MEG_ds.name(3:4), paths.('subj_fcp5'));

% which atlas to use? mmp1 for general, a2009s for aparc fs segmentation
atlas_type = 'a2009s';

% run analysis
fcp_5_restconnectivity(paths, surf_path, atlas_type);
