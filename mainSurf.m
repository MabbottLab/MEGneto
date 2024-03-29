%% Master run file

% Use the following to set up path names and run each step of the pipeline.

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
% startRecordPerformance
analysis_name = 'Right';
project_path = '/mnt/sda/juanita/MEGneto';
rawdata_path = '/mnt/sda/juanita/datasets/right';
mri_path = '/mnt/sda/juanita/MRIs';
paths = megne2setup(project_path, analysis_name, rawdata_path, mri_path, false);
% stopRecordAndDisplay

%% fcp_1: task epoching, jump/muscle artifact detection, bad channel detection
%  To be populated with more information
MEG_ds = struct2table(dir(paths.rawdata));
writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp1'));
<<<<<<< HEAD

=======
>>>>>>> f8ccd90d55d82d3970d91e1f15f13c7eb998f6cc
fcp_1_TaskEpoching(paths)

%% fcp_2: ICA
fcp_2_PreprocessingICA(paths)

%% fcp_2_5: human identifies bad ICA components, reject those components
%  To be populated with more information
fcp_2_5_checkpoint(paths)

%% fcp_3: bad channel repair
%  To be populated with more information
fcp_3_ChannelRepair(paths, included_fcp2)

%%% Import Surfuces from Freesurfer --------------------------
%%% Setup ----------------------------------------------------
%setup paths etc

%sub = 'OIRM02c';
pathsS.basepath = '/data2/DATA_DKI/OIRM_CHASE/meg_freesurfer/'
pathsS.ftpath   = '/home/sonya/matlab/fieldtrip-master'; % this is the path to FieldTrip 

%before running the step below Freesurfer needs to me run and 
% HCP_pipeline (modeifed python script script_ft_postfreesurfer.py) which 
% creates surfaces needed for FieldTrip
fcp_importFreesurferSurfs(pathsS, paths)




%% fcp_4: beamforming
% writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp4'));
%fcp_4_beamforming(paths)

%% fcp_5: connectivity
%fcp_5_taskconnectivity(paths);

