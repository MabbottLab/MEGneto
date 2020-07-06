function trial_summary = getTrialSummary(paths, num_markers, thresh)

fcp1_output = loadjson([paths.anout_grp '/fcp1_output.json']);
ppts = load_participants(paths, 'fcp1');
ppts = ppts.Var1;
ppts = cellfun(@(x) x(1:4), ppts, 'UniformOutput', false);

empty_cells = cellfun(@iscell, fcp1_output.numtrls);
trial_summary(~empty_cells,1) = cell2mat(fcp1_output.numtrls(~empty_cells));
trial_summary(:,2) = num_markers-trial_summary(:,1);
trial_summary(~empty_cells,3) = cell2mat(fcp1_output.HMremove_trls(~empty_cells)); % head motion removed
trial_summary(~empty_cells,4) = cell2mat(fcp1_output.Nremove_trls(~empty_cells)); % noise removed
trial_summary(:,5) = sum(trial_summary(:,2:4),2)./num_markers > thresh;
trial_summary = array2table(trial_summary);
trial_summary.Properties.RowNames = ppts;
trial_summary.Properties.VariableNames = ["num_epochs", ...
                                          "num_incorrect", ...
                                          "hm_removed", ...
                                          "noise_removed", ...
                                          "too_many_removals"];
end