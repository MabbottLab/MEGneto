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

    % set up output variables
    p_val = NaN(size(seed_mat,1),1);
    obs_diffs = p_val;
    eff_size = p_val;
    
    parfor i = 1:size(seed_mat,1) % repeat for each connection between seed and other regions
        if i == this_seed % ignore if it's connectivity with itself
            p_val(i) = NaN;
            obs_diffs(i) = NaN;
            eff_size(i) = NaN;
        else
            % using permutationTest from Github/MATLAB forum
            [p_val(i) obs_diffs(i) eff_size(i)] = permutationTest(seed_mat(i,logical(design_matrix(:,2))), ...
                                seed_mat(i,logical(design_matrix(:,1))), ...
                                num_bootstraps, 'plotresult', 0, 'showprogress', 0);
        end
    end
    % save output by label
    res.(char(atlas.textdata{this_seed,2})) = table(p_val, obs_diffs, eff_size, ...
                                                'VariableNames', ["p_value", "observed_diffs", "effect_size"], ...
                                                'RowNames', string(atlas.textdata(1:90,2)));
end
