function [res, rand_diffs] = bootTestDiffSeeds(paths, seed_regions, freq_band, two_groups, num_bootstraps, thresh, group_names)
% This function performs permutation-based significance testing (via 
% t-test or f-test using the max procedure) to build a null distribution 
% and control for Type 1 error.

% INPUTS-------------------------------------------------------------------
% paths:           - variable containing paths to analysis and data
% seed_regions:    - numeric indices indicating the seed ROIs (e.g. 
%                    if the AAL atlas is used, the default input [1, 2, 3] 
%                    corresponds to the following regions 
%                    ['left precentral gyrus', 'right precentral gyrus', 
%                     'left superior frontal gyrus, dorsolateral']. Note 
%                    that for AAL atlas there are 90 regions, so indices 
%                    should take on values between 1-90). 
% freq_band:       - frequency band of interest
% two_groups:      - true or false to indicate T-max analysis or F-max
%                    analysis, respectively.
% num_bootstraps:  - number of desired bootstraps.
% thresh:          - significance threshold for the p-value.
% group_names:      - array of strings, e.g., ["RAD", "SURG", "TDC"], 
%                     exactly as they appear in the input for group_names 
%                     in the make_NBS function.
% collapsed matrix: - a note on T-max analysis extra input:
% If you wish to perform a T-max analsysis, you must include a .xlsx file
% with two group that will be compared in the T-max analysis for 
% differences. You may use the make_NBS function to prepare this file
% where your input to make_NBS (the ParticipantCategories.xlsx) should have
% two groups you wish to include (e.g. "control" and "surg") and the
% participant IDs that fall under each group. 

% OUTPUTS------------------------------------------------------------------
% res:             - a struct containing the group differences for each
%                    region of interest. p-pos, p-neg and
%                    observed-differences are listed for each. 
% rand_diffs:      - the null distribution 

%% LOAD CONFIG
config      = load_config(paths, paths.name);
config      = config.config;

fprintf("FREQ. BAND: %s\n", freq_band)

%% LOAD ATLAS DATA
if contains(config.beamforming.atlas.filepath, 'mmp') % if MMP glasser atlas
    megneto_path        = fileparts(which('fcp_4_beamforming.m'));
    atlas               = ft_read_atlas([megneto_path '/external/atlas/mmp.mat']);
else
    fullPath                = which('ft_preprocessing.m');
    [pathstr,~,~]           = fileparts(fullPath);
    atlas                   = ft_read_atlas([pathstr config.beamforming.atlas.filepath]);
    if contains(config.beamforming.atlas.filepath, 'aal')
        atlas.tissuelabel   = atlas.tissuelabel(1:90); % we only want non-cerebellar regions (isolate desired regions)
        atlas.tissue(atlas.tissue > 90) = 0;
    end
end

%% CREATE DESIGN MATRIX
load([paths.anout '/NBS' freq_band '_Hz.mat']) % load data corresponding to freq band
load([paths.anout '/design_matrix.mat']) % load design matrix

if two_groups % if a T-max analysis is desired, use collapsed data
    design_matrix = readmatrix([paths.conf_dir '/collapsed_ParticipantCategories.xlsx']);
    rand_diffs = [];
else % if a F-max analysis is desired, use output of make_NBS
    for i = 1:length(design_matrix)
        design_labels(logical(design_matrix(:,i))) = group_names(i);
    end 
    design_labels = cellstr(design_labels);
    res = [];
end

%% ANALYSIS
% for each of the seed regions you requested
for this_seed = seed_regions
    fprintf("ROI: %s\n", atlas.tissuelabel{this_seed});
    seed_mat = squeeze(data_matrix(this_seed, :, :)); % slice the data matrix

    if two_groups % if T-max analysis
        % using permutationTest from Github/MATLAB forum, modified for
        % MaxT procedure
        [p_pos, p_neg, obs_diffs, this_rand_diffs] = permutationTest(seed_mat(:,logical(design_matrix(:,2))), ...
                            seed_mat(:,logical(design_matrix(:,1))), ...
                            num_bootstraps, 'plotresult', 0, 'showprogress', 0);

        % save output by label
        res.(char(atlas.tissuelabel{this_seed})) = table(p_pos', p_neg', obs_diffs', ...
                                                    'VariableNames', ["p_pos", "p_neg", "observed_diffs"], ...
                                                    'RowNames', string(atlas.tissuelabel(setdiff(1:90,this_seed))));
        if any(p_pos < 0.05) || any(p_neg < .05)
            fprintf("\tSignificant connections:\n");
                for i = [find(p_pos < .05) find(p_neg < .05)]
                    fprintf("\t\t%s\n", res.(char(atlas.tissuelabel{this_seed})).Properties.RowNames{i});
                end
        else
            fprintf("\tNo significant connections.\n")
        end

        rand_diffs = [rand_diffs; this_rand_diffs];
    else % if F-max analysis
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
        NOIs{i} = find(strcmp(NOI_origins{i}, atlas.tissuelabel(:)));
        NOIs{i} = [NOIs{i} cell2mat(cellfun(@(x) find(strcmp(x, atlas.tissuelabel(:))), NOI_connections{i}, 'UniformOutput', false))'];
    end
    res.NOIs = NOIs;
end

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
%     res.(char(atlas.tissuelabel{seed_regions(this_seed),2})) = table(p_pos(:,this_seed), ...
%                                                                   p_neg(:,this_seed), ...
%                                                                   obs_diffs(:,this_seed), ...
%                                             'VariableNames', ["p_pos", "p_neg", "observed_diffs"], ...
%                                             'RowNames', string(atlas.tissuelabel(setdiff(1:90,seed_regions(this_seed)),2)));
% end

end
