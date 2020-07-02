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
analysis_name = 'Rest';

project_path = '/home/megneto/data/sbells/meg_freesurfer/MEGneto';

rawdata_path = '/home/megneto/data/sbells/meg_freesurfer/MEG';
mri_path = '/home/megneto/data/sbells/meg_freesurfer/MRIs';
paths = megne2setup(project_path, analysis_name, rawdata_path, mri_path, false);
% stopRecordAndDisplay

%% fcp_1: task epoching, jump/muscle artifact detection, bad channel detection
%  To be populated with more information
MEG_ds = struct2table(dir(paths.rawdata));
writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp1'));

% fcp_1_TaskEpoching(paths)
fcp_1_RestingStateEpoching(paths)

%% fcp_2: ICA
writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp2'));

fcp_2_PreprocessingICA(paths)

%% fcp_2_5: human identifies bad ICA components, reject those components
%  To be populated with more information
writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp2_5'));
fcp_2_5_checkpoint(paths)

%% fcp_3: bad channel repair
%  To be populated with more information
writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp3'));
fcp_3_ChannelRepair(paths)

% fcp_3_ChannelRepair(paths, included_fcp2)


%% fcp_4a: register freesurfer surfaces to ctf mri
% user interaction needed for ac-pc labeling and freesurfer determine
% coordinate (rasa)
writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp4'));
surf_path = '/home/megneto/data/sbells/meg_freesurfer/processed_freesurfer';

fcp_4a_register_CTF_Surf(paths,surf_path)
%no acpac:
%good : 02,03,07,      24,26
%bad  : 06,11,15,16,   17,27,28,
%w acpc (spm12)
%good : 02,03,06,07,11,15,16,17,24,26,27,28,29,30,31,39,41,43,44,47
% 16??17??
%% fcp_4: beamforming surfaces
writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp4'));
% fcp_4_beamforming(paths)

fcp_4_Beamforming_Freq_Surf(paths)
%alpha : (8-13) 10 : 2 
%theta : (4-8) 5 : 2
%beta : 13-30 ; 20 : 7
%gamma (30-48); 39 : 9

%% fcp_5: connectivity
writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp5'));
surf_path = '/home/megneto/data/sbells/meg_freesurfer/processed_freesurfer';
%supported : atlas MMP1.0 (group) : individual a2009s freesufer
atlas_type = 'a2009s'; %'a2009s'mmp1
fcp_5_restconnectivity(paths, surf_path,atlas_type);
