function res = bootTestDiffSeeds(paths, seed_regions, freq_band, num_bootstraps)

atlas = importdata('/mnt/sda/juanita/fieldtrip/template/atlas/aal/ROI_MNI_V4.txt'); % get region labels
load([paths.anout '/NBS_treatmenttype/NBS' freq_band '_Hz.mat']) % load data corresponding to freq band
load([paths.anout '/NBS_treatmenttype/design_matrix.mat']) % load design matrix

% this just collapses the design matrix into patients vs. control since it
% was split into 3
two_groups = true; 
if two_groups
    design_matrix = [sum(design_matrix(:,1:2),2) design_matrix(:,3)];
end

% for each of the seed regions you requested
for this_seed = seed_regions
    seed_mat = squeeze(data_matrix(this_seed, :, :)); % slice the data matrix

    % using permutationTest from Github/MATLAB forum, modified for
    % MaxT procedure
    [p_val, obs_diffs] = permutationTest(seed_mat(:,logical(design_matrix(:,2))), ...
                        seed_mat(:,logical(design_matrix(:,1))), ...
                        num_bootstraps, 'plotresult', 0, 'showprogress', 0);
                    
    % save output by label
    res.(char(atlas.textdata{this_seed,2})) = table(p_val', obs_diffs', ...
                                                'VariableNames', ["p_value", "observed_diffs"], ...
                                                'RowNames', string(atlas.textdata(setdiff(1:90,this_seed),2)));
end
