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
startRecordPerformance
analysis_name = 'TP';
project_path = 'C:\Users\julie tseng\Data';
rawdata_path = 'C:\Users\julie tseng\Data\MEG\left';
mri_path = 'C:\Users\julie tseng\Data\MRIs';
paths = megne2setup(project_path, analysis_name, rawdata_path, mri_path, true);
stopRecordAndDisplay

%% fcp_1: task epoching, jump/muscle artifact detection, bad channel detection
%  To be populated with more information
MEG_ds = struct2table(dir(paths.rawdata));
writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp1'));
fcp_1_TaskEpoching(paths)

%% fcp_2: ICA
%  To be populated with more information
writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp2'));
fcp_2_PreprocessingICA(paths)

%% fcp_2_5: human identifies bad ICA components, reject those components
%  To be populated with more information
fcp_2_5_checkpoint(paths)

%% fcp_3: bad channel repair
%  To be populated with more information
writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp3'));
fcp_3_ChannelRepair(paths)

%% fcp_4: beamforming
% writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp4'));
% fcp_4_beamforming(paths)

% Gutted for pipeline testing
[source, source_proj] = fcp_4_beamforming_gutted(paths);

%% fcp_5: connectivity

% Gutted for pipeline testing
conn_stats = fcp_5_connectivity_gutted(paths, source, source_proj);


