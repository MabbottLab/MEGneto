function [res, rand_diffs] = bootTestDiffSeeds(paths, bootcfg)
% This function performs permutation-based significance testing (via 
% t-test for 2 groups or f-test for >2 groups using the max procedure) 
% to build a null distribution and control for Type 1 error.

% INPUTS-------------------------------------------------------------------
% paths:           - variable containing paths to analysis and data
% bootcfg.
    % seed_regions:    - numeric indices indicating the seed ROIs (e.g. 
    %                    if the AAL atlas is used, the default input [1, 2, 3] 
    %                    corresponds to the following regions 
    %                    ['left precentral gyrus', 'right precentral gyrus', 
    %                     'left superior frontal gyrus, dorsolateral']. Note 
    %                    that for AAL atlas there are 90 regions, so indices 
    %                    should take on values between 1-90). 
    % freq_band:       - frequency band of interest
    % num_bootstraps:  - number of desired bootstraps.
    % thresh:          - significance threshold for the p-value.
    % group_names:      - array of strings, e.g., ["RAD", "SURG", "TDC"], 
    %                     exactly as they appear in the input for group_names 
    %                     in the make_NBS_ready function.
    % nickname:         - nickname of this subanalysis as set in
    %                       make_NBS_ready
%
% OUTPUTS------------------------------------------------------------------
% res:             - a struct containing the group differences for each
%                    region of interest. 
%                    > .p_pos and p_neg represent the proportion of null 
%                      values found above/below the observed statistic
%                    > observed_difference is the actual t- or f-statistic
%                      of the group difference
%                    observed-differences are listed for each. 
% rand_diffs:      - the null distribution 
%
% NOTES--------------------------------------------------------------------
% This function automatically runs a T-max or F-max permutation test based
% on the number of columns (i.e., groups) in the design matrix. 
%
% If you have more than two groups, but want to run a comparison between 
% two groups, you should make a collapsed ParcipantCategories.xlsx sheet
% and re-run make_NBS_ready.m to generate a 2-group design matrix and 
% 2-group data matrices. 

%% LOAD CONFIG
config      = load_config(paths, paths.name);
config      = config.config;

fprintf("FREQ. BAND: %s\n", bootcfg.freq_band)

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
load([paths.anout_grp '/' bootcfg.nickname '_datamat_' bootcfg.freq_band '.mat']) % load data corresponding to freq band
load([paths.anout_grp '/' bootcfg.nickname '_design_matrix.mat']) % load design matrix

if size(design_matrix, 2) == 2 % if only two groups => it's a max T test
    rand_diffs = [];
else % if a F-max analysis is desired, need to generate group name matrix
    for i = 1:length(design_matrix)
        design_labels(logical(design_matrix(:,i))) = bootcfg.group_names(i);
    end 
    design_labels = cellstr(design_labels);
    res = [];
end

%% ANALYSIS
% for each of the seed regions you requested
for this_seed = bootcfg.seed_regions
    fprintf("ROI: %s\n", atlas.tissuelabel{this_seed});
    seed_mat = squeeze(data_matrix(this_seed, :, :)); % slice the data matrix

    if size(design_matrix,2) == 2 % if T-max analysis
        % using permutationTest from Github/MATLAB forum, modified for
        % MaxT procedure
        [p_pos, p_neg, obs_diffs, this_rand_diffs] = permutationTest(seed_mat(:,logical(design_matrix(:,2))), ...
                            seed_mat(:,logical(design_matrix(:,1))), ...
                            bootcfg.num_bootstraps, 'plotresult', 0, 'showprogress', 0);

        % save output by label
        res.(char(atlas.tissuelabel{this_seed})) = table(p_pos', p_neg', obs_diffs', ...
                                                    'VariableNames', ["p_pos", "p_neg", "observed_diffs"], ...
                                                    'RowNames', string(atlas.tissuelabel(setdiff(1:size(seed_mat,1),this_seed))));
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
        res = [res; arrayfun(@(x) anova1(seed_mat(x,:), design_labels, 'off'), 1:(size(seed_mat,1)))];
    end
end

if size(design_matrix,2) == 2 
    NOI_origins = fieldnames(res);
    NOI_connections = cellfun(@(x) ...
                            res.(x).Properties.RowNames(...
                            res.(x){:,'p_pos'} < bootcfg.thresh | res.(x){:,'p_neg'} < bootcfg.thresh), ...
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
