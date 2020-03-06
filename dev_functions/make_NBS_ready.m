function make_NBS_ready(paths, group_names, conn)

% To make data matrix and design matrix for NBS exercise 
% modified from NBS_sam by Sonya - May 9, 2018
% modified by Liz - Sept 2018, Jan 2020
% function-ified by Julie - March 2020

% INPUT--------------------------------------------------------------------
% paths:        the paths struct, as we always input
% group_names:  all group names as an array of strings, e.g., ["surg",
%               "rad", "control"], exactly as they appear in folder names
% conn:         name of connectivity metric as a character array, e.g.,
%               'wpli_deb'
% contrasts:    value that design matrix should take on for each group,
%               thus this should be an array with length (num_groups)

% OUTPUT-------------------------------------------------------------------
% design_matrix.mat:        design matrix for NBS analysis
% NBS[freq]_datamatrix.mat: data matrices organized by frequency band

% IMPORTANT: this version assumes that the participant adjmats only have
% the first entry of the 3rd dimension (the trials dimension) populated.
% This is because WPLI and debiased WPLI must avg across trials for the
% final connectivity estimate. 

%% SETUP PARTICIPANT LISTS

numGroups = length(group_names); % get number of groups
spath = paths.anout; % path to analysis folder

% load info on which participants belong to which groups
groups = readtable([paths.conf_dir '/ParticipantCategories.xlsx']);
groups = groups(:,1:numGroups);

groupPLI = cell(1,numGroups); % for each group
for gg = 1:numGroups 
    this_ppt_list = rmmissing(groups.(group_names(gg))); % isolate ppt list for this group
    groupPLI{gg} = cellfun(@(x) sprintf('%s/%s/fcp_5_adjmat_%s.mat', spath, x, conn), ...
                            this_ppt_list, 'UniformOutput', false); % generate filepaths to those folders
end

%% initialize subjects in conditions - load PLI adj. matrix

nbs_datamat = [];
for gg = 1:numGroups
    fprintf(['\t---- Group ',num2str(gg),' ----\n']);
        for ss = 1:length(groupPLI{gg})
            this_adjmat = load(groupPLI{gg}{ss}); % load the struct
            nbs_datamat = cat(3, nbs_datamat, this_adjmat.adjmat(:,:,1,:)); % throw in the participant's conn
        end
end

%% Save matrix by frequency band 

freq ={'theta','alpha','beta','Lgamma','Hgamma'};
num_freqs = size(nbs_datamat,4); % just in case not every freq is there

for ff = 1:num_freqs % range freq save
    data_matrix = nbs_datamat(:,:,:,ff);
    datamat_output = [spath,'/NBS',freq{ff},'_Hz.mat'];
    save(datamat_output,'data_matrix','-mat');
    fprintf(['NBS Data Matrix saved to: ',datamat_output,'\n']);
end

%% make design matrix for NBS
% modified from Julie Sato's NBS script (written by Simeon Wong)

design_matrix = [];
for gg = 1:numGroups
    this_gp = zeros(1,numGroups);
    this_gp(gg) = 1;
    design_matrix = vertcat(design_matrix, repmat(this_gp, length(groupPLI{gg}), 1));
end

% save 
designmat_output = [spath,'/design_matrix.mat'];
save(designmat_output,'design_matrix','-mat');
fprintf(['\nDesign Matrix saved to: ',designmat_output,'\n\n']);

end
