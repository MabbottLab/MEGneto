function [out, failure] = freeviewing_ds_pid_match(paths,step)
% Compared to the standard MEG pipeline, the following lines were edited:
% 20-30

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

% if length(pid) < length(out.ds) % check for mismatched MEG/MRI data
%     fprintf('\n\n%s\n', 'Subjects missing MRI files:')
%     TF = contains(out.ds, pid);
%     for i = 1:length(TF)
%         if TF(i) == 0
%             disp(out.ds(i))
%         end
%     end
%     error('The participants above are missing MRI data. Please investigate, edit subj_%s.csv, and re-run %s.', ...
%          step, step)
% end

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