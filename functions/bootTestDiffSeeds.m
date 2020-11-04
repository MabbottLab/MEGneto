function [res, rand_diffs] = bootTestDiffSeeds(paths, seed_regions, freq_band, two_groups, num_bootstraps, thresh)

fprintf("FREQ. BAND: %s\n", freq_band)

atlas = importdata('/mnt/sda/juanita/fieldtrip/template/atlas/aal/ROI_MNI_V4.txt'); % get region labels
load([paths.anout '/NBS_treatmenttype/NBS' freq_band '_Hz.mat']) % load data corresponding to freq band
load([paths.anout '/NBS_treatmenttype/design_matrix.mat']) % load design matrix

% this just collapses the design matrix into patients vs. control since it
% was split into 3
if two_groups
    design_matrix = [sum(design_matrix(:,1:2),2) design_matrix(:,3)];
    rand_diffs = [];
else
    design_labels(logical(design_matrix(:,1))) = "rad";
    design_labels(logical(design_matrix(:,2))) = "surg";
    design_labels(logical(design_matrix(:,3))) = "tdc";
    design_labels = cellstr(design_labels);
    res = [];
end

% for each of the seed regions you requested
for this_seed = seed_regions
    fprintf("ROI: %s\n", atlas.textdata{this_seed,2});
    seed_mat = squeeze(data_matrix(this_seed, :, :)); % slice the data matrix

    if two_groups
        % using permutationTest from Github/MATLAB forum, modified for
        % MaxT procedure
        [p_pos, p_neg, obs_diffs, this_rand_diffs] = permutationTest(seed_mat(:,logical(design_matrix(:,2))), ...
                            seed_mat(:,logical(design_matrix(:,1))), ...
                            num_bootstraps, 'plotresult', 0, 'showprogress', 0);

        % save output by label
        res.(char(atlas.textdata{this_seed,2})) = table(p_pos', p_neg', obs_diffs', ...
                                                    'VariableNames', ["p_pos", "p_neg", "observed_diffs"], ...
                                                    'RowNames', string(atlas.textdata(setdiff(1:90,this_seed),2)));
        if any(p_pos < 0.05) || any(p_neg < .05)
            fprintf("\tSignificant connections:\n");
                for i = [find(p_pos < .05) find(p_neg < .05)]
                    fprintf("\t\t%s\n", res.(char(atlas.textdata{this_seed,2})).Properties.RowNames{i});
                end
        else
            fprintf("\tNo significant connections.\n")
        end

        rand_diffs = [rand_diffs; this_rand_diffs];
    else
        % end up with: (num_seeds) x (paired_region) x (frequency band) matrix
        res = [res; arrayfun(@(x) anova1(seed_mat(x,:), design_labels, 'off'), 1:90)];
    end
end

if two_groups
    NOI_origins = fieldnames(res);
    NOI_connections = cellfun(@(x) ...
                            res.(x).Properties.RowNames(...
                            res.(x){:,'p_pos'} < thresh | res.(x){:,'p_neg'} < thresh), ...
                            NOI_origins, 'UniformOutput', false);
    for i = 1:length(NOI_origins)
        NOIs{i} = find(strcmp(NOI_origins{i}, atlas.textdata(:,2)));
        NOIs{i} = [NOIs{i} cell2mat(cellfun(@(x) find(strcmp(x, atlas.textdata(:,2))), NOI_connections{i}, 'UniformOutput', false))'];
    end
end
res.NOIs = NOIs;

%%% COLLAPSING ACROSS ALL SEED REGIONS

% num_ppts = size(data_matrix,3);
% num_seeds = length(seed_regions);
% num_conns = size(data_matrix,1);
% seed_mat = reshape(permute(squeeze(data_matrix(seed_regions,:,:)), [2 1 3]), num_seeds*num_conns, num_ppts);
% 
% % using permutationTest from Github/MATLAB forum, modified for
% % MaxT procedure
% [p_pos, p_neg, obs_diffs] = permutationTest(seed_mat(:,logical(design_matrix(:,2))), ...
%                     seed_mat(:,logical(design_matrix(:,1))), ...
%                     num_bootstraps, 'plotresult', 0, 'showprogress', 0);
%                 
% p_pos = reshape(p_pos, num_conns-1, num_seeds);
% p_neg = reshape(p_neg, num_conns-1, num_seeds);
% obs_diffs = reshape(obs_diffs, num_conns-1, num_seeds);
%                     
% % save output by label
% for this_seed = 1:num_seeds
%     res.(char(atlas.textdata{seed_regions(this_seed),2})) = table(p_pos(:,this_seed), ...
%                                                                   p_neg(:,this_seed), ...
%                                                                   obs_diffs(:,this_seed), ...
%                                             'VariableNames', ["p_pos", "p_neg", "observed_diffs"], ...
%                                             'RowNames', string(atlas.textdata(setdiff(1:90,seed_regions(this_seed)),2)));
% end

end
