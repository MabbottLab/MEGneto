%% Master run file

% Use the following to set up path names and run each step of the pipeline.

%% if you need to re-load paths
addpath(genpath('/mnt/sda/juanita/MEGneto'))%% if you need to re-load paths
addpath(genpath('/mnt/sda/juanita/MEGneto'))
addpath('/mnt/sda/juanita/fieldtrip')
ft_defaults
paths = loadjson('mnt/sda/juanita/MEGneto/analysis/left/config/paths.json');
names = fieldnames(paths);
for i = 1:100
    paths.(char(names{i})) = ['/mnt/sda/juanita/MEGneto/' paths.(char(names{i}))];
end

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
analysis_name = 'left';
project_path = '/mnt/sda/juanita/MEGneto';
rawdata_path = '/mnt/sda/juanita/datasets/left';
mri_path = '/mnt/sda/juanita/MRIs';
paths = megne2setup(project_path, analysis_name, rawdata_path, mri_path, false);
% stopRecordAndDisplay

%% fcp_1: task epoching, jump/muscle artifact detection, bad channel detection
%  To be populated with more information
MEG_ds = struct2table(dir(paths.rawdata));
writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp1'));
fcp_1_TaskEpoching(paths)

%% fcp_2: ICA
fcp_2_PreprocessingICA(paths)

%% fcp_2_5: human identifies bad ICA components, reject those components
%  To be populated with more information
fcp_2_5_checkpoint(paths)

%% fcp_3: bad channel repair
%  To be populated with more information
fcp_3_ChannelRepair(paths, included_fcp2)

%% fcp_4: beamforming
% writecell(MEG_ds.name(3:(height(MEG_ds))), paths.('subj_fcp4'));
fcp_4_beamforming(paths)

%% fcp_5: connectivity
fcp_5_taskconnectivity(paths);

%% prepare for NBS

group_names = ["RAD","SURG", "TDC"];
conn = 'wpli_debiased';
make_NBS_ready(paths, group_names, conn)

%% prepare for bnv

addpath('/home/liz/matlab/BCT/2016_01_16_BCT');
addpath('/home/liz/matlab/BrainNetViewer')

brainnet = [];
brainnet.netmetric = 'nbs'; % NETWORK METRIC ('pls' or 'nbs')
brainnet.nbs_type = 'extent';
brainnet.connmetric = 'pli'; % connectivity metric
brainnet.fb_interest = 'Lgamma'; % freq band
brainnet.file_name = 'NBS_RAD_vs_SURG'; % name of output
brainnet.edge_weight = 1;
brainnet.colour = 2; 
brainnet.size_nodes = 1;
brainnet.grp = {1:19, 20:32};
brainnet.noi = [19,29,36,68,83,85];

make_BNV_ready(paths, brainnet)



