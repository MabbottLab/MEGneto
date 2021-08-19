function paths = freeviewing_megne2setup(project_path, analysis_name, rawdata_path, mri_path, overwrite)

% Compared to the standard MEG pipeline, the following lines were edited:
% 55

% MEGNE2SETUP will create a subfolder to PROJECT_PATH or the current 
% working directory (if PROJECT_PATH is not provided) called ANALYSIS_NAME, 
% create config and analysis subfolders within, and create unfilled 
% participants.txt and config.json files in the config directory. If you 
% already have a config file with your preferred parameters, replace the 
% unfilled config file with that one. You should have a rawdata folder in 
% your project directory with all the .ds files, and a subfolder to rawdata
% called "MRIs" which contains all the MRIs in the form [PID].mri.
%
% NOTES:
%   The top-level fieldtrip and MEGneto folders should be visible on the
%   path at the outset of this code. Subfolders will be added accordingly. 
% 
% INPUTS:
%   project_path        =   string: name and location of project directory 
%                           (if no location defined, then it will be created 
%                           in current working directory. 
%   analysis_name       =   string
%   rawdata_path        =   string defining MEG data path (*.ds folders)
%   mri_path            =   string defining MRI data path (*.mri)
%   overwrite           =   boolean t/f (default = false)
%
% RETURNS:
%   paths               =   struct defining paths to data, participant
%                           folders, analysis folders, config files, etc. 
%
% See also: PATH_GENERATION, PATH_CHECK

% Last updated by: Julie Tseng, 2020-01-07
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.

%% SET UP LOGGING FILE

right_now = clock;
mkdir([project_path '/analysis/' analysis_name '/config/'])
log_filename = [project_path '/analysis/' analysis_name '/config/log_' ...
    sprintf('%d%d%d', right_now(1:3))];
diary(log_filename)

fprintf('%d:%d:%02.f       Now running **%s**.\n', ...
    right_now(4:6), mfilename)

%% GET IMPORTANT PATHS AND ADD APPROPRIATE FOLDERS

% identify where MEGneto is running
megne2_name = which(mfilename());
megne2_name = megne2_name(1:end-(length(['/' mfilename() '.m'])));

megne2_name = '/mnt/hpc-megneto/MEGneto';

prior_dir = pwd();      % Keep track of prior directory
cd(megne2_name)         % Set path to main MEGne code
addpath('.')
addpath(genpath('functions'))   % Add path of functions
addpath('external/glob')        % glob makes finding paths very easy
addpath('external/jsonlab')     % json parser

% check whether fieldtrip is visible
addpath(strcat(prior_dir, '/fieldtrip'))
if ~exist('ft_defaults','file')                 % If not, throw error
    error('FieldTrip not visible in path.')
else
    ft_defaults;                                % If so, run FT setup file
end

cd(prior_dir);

%% CHECK INPUT PARAMETERS

% set the overwrite default value if not specified
if nargin < 5
    overwrite = false;
end

% check if project path is empty or non-existent
if strcmp(project_path,'') || ~exist('project_path', 'var') 
    if ~strcmp(megne2_name, pwd())
        project_path = pwd();
        warning('No project_path provided. Default set to current working directory.\n')
    else
        error('No project_path provided and current working directory is MEGne code directory.')
    end
end

% standardizes path to end without /
if strcmp(project_path(end),'/') 
    project_path = project_path(1:end-1);
end

% checks for existence, fullness, and correct typing of the analysis name
if ~ischar(analysis_name) || strcmp(analysis_name,'') 
    error('Analysis name must be non-empty value of class char')
end

% checks for existence of MEG data path
if ~exist(rawdata_path,'dir')
    error('Rawdata directory does not exist')
end

% checks for existence of MRI data path
if ~exist(mri_path,'dir')
    error('MRIs directory does not exist')
end

% checks whether there already exists an analysis output folder and no
% indication that it should be overwritten
if exist([project_path '/' analysis_name], 'file') && overwrite == false
    error(['A directory or file of name "' analysis_name '" already exists. Rename it, remove it or change the name of this analysis before trying again.'])
end

%% MAKE DIRECTORIES AND PATHS STRUCT

% run path_generation to generate locations
[paths, subj] = freeviewing_path_generation(project_path, analysis_name, rawdata_path, mri_path, overwrite);

% FOLDER STRUCTURE
% project_path
%   -- analysis
%       -- analysis_name #1
%           -- analysis
%               -- group
%               -- ST01
%               -- ...
%               -- STXX
%           -- config
%       -- analysis_name #2
%           -- analysis
%               -- group
%               -- ST01
%               -- ...
%               -- STXX
%           -- config

% create top-level analysis output folder that holds multiple analyses
if ~exist(paths.analyses, 'dir')
    mkdir(paths.analyses);
end
mkdir(paths.anhome);            % create analysis-specific home folder
mkdir(paths.conf_dir);          % create analysis-specific config folder
mkdir(paths.anout);             % create analysis-specific analysis folder
mkdir(paths.anout_grp);         % create analysis-specific group folder

% initialize JSON config file based on existence, template, overwrite
if overwrite == true || ~exist(paths.mainconf,'file')
    if exist([megne2_name '/templates/' analysis_name '.json'],'file')
        copyfile([megne2_name '/templates/' analysis_name '.json'],...
            [paths.conf_dir '/' analysis_name '.json'])
    else
        copyfile([megne2_name '/configs/empty_config.json'],...
            [paths.conf_dir '/' analysis_name '.json'])
    end
end

% loop over number of subjects to make analysis directories for each
for ii = 1:height(subj)
    mkdir(char(paths.(subj.pids{ii})));
end

% create fcp csvs that facilitate participant inclusion/exclusion
% these files are filled with PIDs by the user prior to each fcp step
paths_cell = table2cell(paths)';
csv = ~cellfun(@isempty, (strfind(paths_cell,'csv')));
csv = paths_cell(csv);
for ii = 1:length(csv)
    if ~exist(csv{ii}, 'file')
        fid = fopen(csv{ii},'wt');
        fclose(fid);
    end
end
writetable(subj, paths.all_subj_pids);

% save JSON paths for later reloading
save_to_json(paths,paths.paths); 

%%  run final check that everything initialized correctly
path_check(paths, {'name'});

%% turn off logging

right_now = clock;
fprintf('%d:%d:%02.f       Done running **%s**.\n', ...
    right_now(4:6), mfilename)
diary off

