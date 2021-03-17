% Modify a paths struct and re-save for moving to a different computer
% assuming that you are starting from fcp_4_beamforming

% change project path locations
new_project_path = '/home/jtseng/katie';
paths_loc = '/home/jtseng/katie/analysis/SpeedofThinkingMEG/config/paths.json';
paths = loadjson(paths_loc);
fn = fieldnames(paths);
to_replace = paths.(fn{1});

for k = 1:numel(fn)
    paths.(fn{k}) = strrep(paths.(fn{k}), to_replace, new_project_path);
end

% change raw MRI path
paths.rawmri = '/home/jtseng/noor/MEG_T1s';

% save new json
save_to_json(paths, paths.paths);