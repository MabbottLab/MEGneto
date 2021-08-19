function freeviewing_write_match_if_not_empty(paths,step)

% Compared to the standard MEG pipeline, the following lines were edited:
% 11

%WRITE_MATCH_IF_NOT_EMPTY accepts a paths structure and step argument of
%class char and, if the provided table contains empty values, backs up the
%previous before saving the new one.

subj_ds = load_participants(paths,step);
[subj_match, failure] = freeviewing_ds_pid_match(paths,step);
if isempty(subj_match) || failure
    error('No participants selected / participants selected improperly')
end
if check_csv_has_empty(paths.(['subj_' step '_match']));
    warning(['current saved' step ' ds-pid match table is empty'...
        'or has unmatched values, backing up and overwriting. '...
        'This is normal on a first run.'])
    copyfile(paths.(['subj_' step '_match']), [paths.(['subj_' step '_match']) '.bak'])
end
writetable(subj_match,paths.(['subj_' step '_match']));