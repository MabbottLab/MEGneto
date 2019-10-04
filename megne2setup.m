function paths = megne2setup(project_path, analysis_name, rawdata_path, mri_path, overwrite)
%%megne2setup
% This function will create a subfolder to PROJECT_PATH
% or the current working directory if PROJECT_PATH is not provided called
% ANALYSIS_NAME, create config
% and analysis subfolders to that, and create unfilled participants.txt and
% config.json files in the config directory. If you already have a config
% file with your preferred parameters, replace the unfilled config file with
% that one. You should have a rawdata folder in your project directory with
% all the .ds files, and a subfolder to rawdata called "MRIs" which
% contains all the MRIs in the form [PID].mri.

%% INPUT

% Identify where megne2 is running
megne2_name = which(mfilename());
megne2_name = megne2_name(1:end-(length(['/' mfilename() '.m'])));
if nargin < 5
    overwrite = false;
end

%% PATH SETUP

prior_dir = pwd();
cd(megne2_name)
addpath('.')
addpath(genpath('functions'))
addpath('external/glob')
addpath('external/jsonlab')
fieldtrip = glob('external/fieldtrip*');
if length(fieldtrip) > 1; warning(['Multiple fieldtrip versions, using ' fieldtrip{1}]); end
addpath(fieldtrip{1});
ft_defaults();
cd(prior_dir);

%% INPUT CHECKING AND CLEANING
if strcmp(project_path,'') || ~exist('project_path', 'var') %Checks if project path is empty or non-existent
    if ~strcmp(megne2_name, pwd())
        project_path = pwd();
        warning('No project_path provided. Defaulting to current working directory')
    else
        error('No project_path provided and current working directory is MEGne2 code directory')
    end
end
if strcmp(project_path(end),'/') % Standardizes path to end without /
    project_path = project_path(1:end-1);
end
if ~ischar(analysis_name) || strcmp(analysis_name,'') % Checks for existence, fullness, and correct typing of the analysis name
    error('NAME must be non-empty value of class char')
end
if ~exist(rawdata_path,'dir')
    error('Rawdata directory does not exist')
end
if ~exist(mri_path,'dir')
    error('MRIs directory does not exist')
end
if exist([project_path '/' analysis_name], 'file') && overwrite == false
    error(['A directory or file of name "' analysis_name '" already exists. Rename it, remove it or change the name of this analysis before trying again.'])
end

%% PATH SETUP

[paths, subj] = path_generation(project_path, analysis_name, rawdata_path, mri_path);

%% FILE CREATION

if ~exist(paths.analyses, 'dir')
    mkdir(paths.analyses);
end
mkdir(paths.anhome);
mkdir(paths.conf_dir);
mkdir(paths.anout);
mkdir(paths.anout_grp);

save_to_json(paths,paths.paths);
if overwrite == true || ~exist(paths.mainconf,'file')
    if exist([megne2_name '/templates/' analysis_name '.json'],'file')
        copyfile([megne2_name '/templates/' analysis_name '.json'],...
            [paths.conf_dir '/' analysis_name '.json'])
    else
        copyfile([megne2_name '/configs/empty_config.json'],...
            [paths.conf_dir '/' analysis_name '.json'])
    end
end

for ii = 1:height(subj)
    mkdir(paths.(subj.pids{ii}));
end

paths_cell = table2cell(paths)';
csv = ~cellfun(@isempty, (strfind(paths_cell,'csv')));
csv = paths_cell(csv);
for ii = 1:length(csv)
    if ~exist(csv{ii}, 'file')
        fid = fopen(csv{ii},'wt');
        fclose(fid);
    end
end

writetable(subj,paths.all_subj_pids);
path_check(paths, {'name'});

