function [paths, pids] = path_generation(project_path, analysis_name, rawdata_path, mri_path)

% PATH_GENERATION sets up the locations of all relevant input files in a 
% depth-1 table. This struct feeds forward into each fcp_x step. If a paths
% struct has already been saved as a JSON to the appropriate location, this
% code will retrieve that struct and return it. 
% 
% INPUTS:
%   project_path        =   string: name and location of project directory 
%                           (if no location defined, then it will be created 
%                           in current working directory. 
%   analysis_name       =   string
%   rawdata_path        =   string defining MEG data path (*.ds folders)
%   mri_path            =   string defining MRI data path (*.mri)
%
% RETURNS:
%   paths               =   struct defining paths to data, participant
%                           folders, analysis folders, config files, etc. 
%   pids                =   table of all participant IDs
%
% See also: MEGNE2SETUP, PATH_CHECK

% Last updated by: Julie Tseng, 2020-01-07
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.
%   

%% Begin code

paths.home = project_path;      % home folder of project 
paths.analyses = [paths.home '/analysis'];              % umbrella analysis folder
paths.anhome = [paths.analyses '/' analysis_name]; % specific analysis folder
paths.conf_dir = [paths.anhome '/config'];  % specific analysis config files
paths.paths = [paths.conf_dir '/paths.json']; % location of JSON record of paths struct

% if this paths struct has already been generated before
if exist(paths.paths, 'file') % then retrieve it and return
    warning(['paths.json already exists, retrieving from ' paths.paths]);
    paths = loadjson(paths.paths);
    paths = struct2table(paths);
    pids = readtable(paths.all_subj_pids);
else % otherwise, define the rest of the paths struct
    paths.name = analysis_name;     % analysis name (e.g., BOTH)
    paths.anout = [paths.anhome '/analysis'];               % specific analysis OUTPUT folder
    paths.anout_grp = [paths.anout '/group'];               % specific analysis GROUP output folder
    paths.mainconf = [paths.conf_dir '/' paths.name '.json'];       % location of JSON config file
    paths.all_subj_pids = [paths.conf_dir '/all_subj_pids.csv'];    % CSV of ALL participants under analysis

    step_nums = ["1", "2", "2_5", "3", "4", "5"];
    for val = step_nums
        % initialize empty CSVs to define participants included at each stage
        % to be filled manually by the user
        paths.(sprintf('subj_fcp%s',val)) = [paths.conf_dir '/subj_fcp' char(val) '.csv'];
        
        % initialize empty CSVs to define *matched* participants who have both MRI
        % and MEG datasets        
        paths.(sprintf('subj_fcp%s_match',val)) = [paths.conf_dir '/subj_match_fcp', char(val), '.csv'];
    end

    % paths to MEG and MRI data
    paths.rawdata = rawdata_path;
    paths.rawmri = mri_path;

    % individual MRI data paths
    paths = struct2table(paths);
    all_participants = glob([paths.rawmri '/*.mri']);
    pids = cell(length(all_participants),1);
    for ii = 1:length(all_participants)
        if contains(all_participants{ii}, '_') % check whether the character array contains an underscore 
            extracted_pid = extractBetween(all_participants{ii}, 'MRIs/', '_'); % extract characters between 'MRIs/' and '_'
        else
            extracted_pid = extractBetween(all_participants{ii}, 'MRIs/', '.'); % extract characters between 'MRIs/' and '.'
        end 
        pids{ii} = extracted_pid{1}; % extract and store content of 1x1 cell array containing the pid 
        all_participants{ii} = [paths.anout '/' extracted_pid{1}];
        
    end
    all_participants = cell2struct(all_participants,pids);
    all_participants = struct2table(all_participants);

    % append individual participant paths to overall paths struct
    paths = [paths, all_participants]; 

    % save and return all PIDs
    pids = cellfun(@char,pids,'UniformOutput',false); 
    pids = cell2table(pids);
end

end

