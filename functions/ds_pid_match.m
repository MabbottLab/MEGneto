function [out, failure] = ds_pid_match(paths,step)
%DS_PID_MATCH extracts the list of all participants, matches the list of
%PIDs from the MRI files, and outputs a csv with each input DS matched to a
%participant ID. In cases where there is ambiguity, it asks you to review
%and correct the matching and returns failure = true
pid = readtable(paths.all_subj_pids);
pid = pid.pids;
ds = load_participants(paths, step);
for ii = 1:length(pid)
    pid{ii} = char(pid{ii});
end
out.ds = ds.Var1;
out.pid = cell(height(ds),1);
failure = false;
for ii = 1:height(ds)
    match = cellfun(@(pid) strfind(out.ds{ii},pid), pid, 'UniformOutput', false);
    match = cellfun(@(x) ~isempty(x),match);
    match = find(match);
    if length(match) == 1
        out.pid{ii} = pid(match);
    elseif isempty(match)
        warning(['Automatic matching was unable to match ' ds.Var1{ii}...
            ' to a subject id (None found). Please fill it out manually '...
            'in the ' step ' match csv file in the config']);
        failure = true;
    elseif length(match) > 1
        warning(['Automatic matching was unable to match ' ds.Var1{ii}...
            ' to a subject id (Multiple found). Please fill it out manually '...
            'in the ' step ' match csv file in the config']);
        failure = true;
    end
end
out = struct2table(out);
out.pid = cellfun(@char,out.pid,'UniformOutput',false);