function [paths, pids] = path_generation(project_path, analysis_name, rawdata_path, mri_path)
%PATH_GENERATION Path struct generator for Megne2 pipeline
%   Sets up the locations of all relevant input files in a depth-1 table.
paths.name = analysis_name;

paths.home = project_path;
paths.analyses = [paths.home '/analysis'];
paths.anhome = [paths.analyses '/' analysis_name];
paths.anout = [paths.anhome '/analysis'];
paths.anout_grp = [paths.anout '/group'];
paths.conf_dir = [paths.anhome '/config'];
paths.mainconf = [paths.conf_dir '/' paths.name '.json'];
paths.paths = [paths.conf_dir '/paths.json'];
paths.all_subj_pids = [paths.conf_dir '/all_subj_pids.csv'];
paths.subj_fcp1 = [paths.conf_dir '/subj_fcp1.csv'];
paths.subj_fcp2 = [paths.conf_dir '/subj_fcp2.csv'];
paths.subj_fcp3 = [paths.conf_dir '/subj_fcp3.csv'];
paths.subj_fcp4 = [paths.conf_dir '/subj_fcp4.csv'];
paths.subj_fcp5 = [paths.conf_dir '/subj_fcp5.csv'];
paths.subj_fcp1_match = [paths.conf_dir '/subj_match_fcp1.csv'];
paths.subj_fcp2_match = [paths.conf_dir '/subj_match_fcp2.csv'];
paths.subj_fcp3_match = [paths.conf_dir '/subj_match_fcp3.csv'];
paths.subj_fcp4_match = [paths.conf_dir '/subj_match_fcp4.csv'];
paths.subj_fcp5_match = [paths.conf_dir '/subj_match_fcp5.csv'];
paths.rawdata = rawdata_path;
paths.rawmri = mri_path;

paths = struct2table(paths);

all_participants = glob([paths.rawmri '/*.mri']);
pids = cell(length(all_participants),1);
for ii = 1:length(all_participants)
    pids{ii} = all_participants{ii}(length(paths.rawmri)+2:end-7);
    all_participants{ii} = [paths.anout '/' all_participants{ii}(length(paths.rawmri)+2:end-7)];
end
all_participants = cell2table(all_participants', 'VariableNames', pids);

paths = [paths, all_participants];
pids = cellfun(@char,pids,'UniformOutput',false);
pids = cell2table(pids);

if exist(paths.paths, 'file')
    warning(['paths.json already exists, retrieving from ' paths.paths]);
    paths = loadjson(paths.paths);
    paths = struct2table(paths);
end
end

