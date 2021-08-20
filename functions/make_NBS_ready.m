function make_NBS_ready(paths, group_names, conn, nickname)

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
% nickname:     nickname of subanalysis prepended to function outputs

% OUTPUT-------------------------------------------------------------------
% [nickname]_designmatrix.mat:          design matrix for NBS analysis
% [nickname]_datamatrix_[freq].mat:     data matrices organized by frequency band.
%                                       Note that the order of participants in these
%                                       output files is re-arranged - it matches the
%                                       order of participants in the design matrix
%                                       which is clarified in the output below
%                                       (ordered_ppt_list.txt).
% [nickname]_orderedpptlist.txt:        list of participants in order of the design
%                                       matrix rows for reference purposes

% IMPORTANT: this version assumes that the participant conn_mats only have
% the first entry of the 3rd dimension (the trials dimension) populated.
% This is because WPLI and debiased WPLI must avg across trials for the
% final connectivity estimate. 

%% LOAD JSON CONFIG FILE
config      = load_config(paths, paths.name);
config      = config.config;

%% SETUP PARTICIPANT LISTS

numGroups = length(group_names); % get number of groups
spath = paths.anout; % path to analysis folder
spath_grp = paths.anout_grp;

% load info on which participants belong to which groups
% groups = readtable([paths.conf_dir '/' nickname '_ParticipantCategories.xlsx']);
% groups = groups(:,1:numGroups);
% groupPLI = cell(1,numGroups); % for each group
% ppt_list = cell(1,numGroups); % to keep track of all ppts
% for gg = 1:numGroups
%     this_ppt_list = rmmissing(groups.(group_names(gg))); % isolate ppt list for this group
%     groupPLI{gg} = cellfun(@(x) sprintf('%s/%s/fcp_5_conn_mat_%s.mat', spath, x, conn), ...
%                        config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp5';

subj_match  = ds_pid_match(paths,step);
ssSubjPath  = @(x) paths.(subj_match.pid{x});

for trial = 24
    catmatrix_out = [];
    for ss = 1:length(subj_match.ds)
        try
            load([ssSubjPath(ss) '/NANcompiledCATmatrix.mat'], '-mat');
        catch
            load([ssSubjPath(ss) '/AAL_beamforming_results.mat'], '-mat');
        end
        if length(newcatmatrix) < length(catmatrix_out)
            newcatmatrix(length(newcatmatrix)+1:length(catmatrix_out), :, :) = NaN;
        elseif length(newcatmatrix) > length(catmatrix_out)
            catmatrix_out(length(catmatrix_out)+1:length(newcatmatrix), :, :) = NaN;
        end
        if ss == 1
            catmatrix_out = newcatmatrix(:, trial, :);
        else
            catmatrix_out = cat(2, catmatrix_out, newcatmatrix(:, trial, :));
        end
    end
    save([paths.anout '/trial/trial' num2str(trial) 'compiledCATmatrix' '.mat'],'catmatrix_out','-v7.3');
end
     this_ppt_list, 'UniformOutput', false); % generate filepaths to those folders
%     ppt_list{gg} = this_ppt_list; % add current group's ppt list to over ppt list
% end
% 
% ordered_ppt_list = [];
% for i = 1:length(ppt_list)
%     ordered_ppt_list = cat(1, ordered_ppt_list, cell2mat(ppt_list{i}));
% end
% 
% % save ordered participant list
% orderedpptlist_output = [spath_grp,'/' nickname '_orderedpptlist.txt'];
% dlmwrite(orderedpptlist_output, ordered_ppt_list, '');
% fprintf(['\nParticipant list saved to: ',orderedpptlist_output,'\n\n']);

%% initialize subjects in conditions - load PLI conn. matrix
step        = 'fcp5';
subj_match  = ds_pid_match(paths,step);
ssSubjPath  = @(x) paths.(subj_match.pid{x});

groupPLI = cell(1, numGroups);
for ss = 1:length(subj_match.ds)

    for gg = 1:numGroups
        if gg == 1
            groupPLI{gg} = {[ssSubjPath(ss) '/fcp_5_controls_conn_mat_wpli_debiased.mat'], '-mat'}; % generate filepaths to those folders
        else
            groupPLI{gg} config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp5';

subj_match  = ds_pid_match(paths,step);
ssSubjPath  = @(x) paths.(subj_match.pid{x});

for trial = 24
    catmatrix_out = [];
    for ss = 1:length(subj_match.ds)
        try
            load([ssSubjPath(ss) '/NANcompiledCATmatrix.mat'], '-mat');
        catch
            load([ssSubjPath(ss) '/AAL_beamforming_results.mat'], '-mat');
        end
        if length(newcatmatrix) < length(catmatrix_out)
            newcatmatrix(length(newcatmatrix)+1:length(catmatrix_out), :, :) = NaN;
        elseif length(newcatmatrix) > length(catmatrix_out)
            catmatrix_out(length(catmatrix_out)+1:length(newcatmatrix), :, :) = NaN;
        end
        if ss == 1
            catmatrix_out = newcatmatrix(:, trial, :);
        else
            catmatrix_out = cat(2, catmatrix_out, newcatmatrix(:, trial, :));
        end
    end
    save([paths.anout '/trial/trial' num2str(trial) 'compiledCATmatrix' '.mat'],'catmatrix_out','-v7.3');
end
= cellfun(@(x) sprintf('%s/%s/fcp_5_%s_conn_mat_%s.mat', spath, x, 'social', conn), ...
                'UniformOutput', false); % generate filepaths to those folders
        end
    end
end

nbs_datamat = [];

for gg = 1:numGroups
    fprintf(['\t---- Group ',num2str(gg),' ----\n']);
        for ss = 1:length(groupPLI{gg})
            this_conn_mat = load(groupPLI{gg}{ss}); % load the struct
            this_conn_mat = this_conn_mat.conn_mat;
            nbs_datamat = cat(4, nbs_datamat, this_conn_mat); % throw in the participant's conn
        end
end

%% Save matrix by frequency band 

freq = config.connectivity.freq_names;
nbs_datamat = permute(nbs_datamat, [1 2 4 3]);
num_freqs = size(nbs_datamat,4); % just in case not every freq is there

for ff = 1:num_freqs % range freq save
    data_matrix = nbs_datamat(:,:,:,ff);
    datamat_output = [spath_grp,'/' nickname '_datamatrix_' freq{ff} '.mat'];
    save(datamat_output,'data_matrix','-mat');
    fprintf(['NBS Data Matrix saved to: ',datamat_output,'\n']);
end

%% make design matrix for NBS

design_matrix = [];
for gg = 1:numGroups
    this_gp = zeros(1,numGroups);
    this_gp(gg) = 1;
    design_matrix = vertcat(design_matrix, repmat(this_gp, length(groupPLI{gg}), 1));
end

%% make design matrix for dorsal prefrontal cortex, inferior parietal lobe, FEF, and PEF
indices = [18, 19, 20, 21, 10, 29]';
tmp = indices + 50;
indices = [indices;tmp];
data_matrix = cat(3, social(indices, indices, :, 2), control(indices, indices, :, 2));
data_matrix = nbs_datamat(:, :, :, 1);
% save 
designmat_output = [spath_grp,'/' nickname '_designmatrix.mat'];
save(designmat_output,'design_matrix','-mat');
fprintf(['\nDesign Matrix saved to: ',designmat_output,'\n\n']);

end
