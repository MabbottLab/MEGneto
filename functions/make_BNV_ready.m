function make_BNV_ready(paths, brainnet)

% To make *.node and *.edge files for viewing connectivity results from NBS 
% or custom max-T analysis on BrainNet Viewer (BNV). 
%
% PLEASE NAME the saved NBS workspace variable as
% [nickname]_NBSResults_[freq].mat
% 
% INPUT--------------------------------------------------------------------
% paths:        the paths struct, as we always input
% brainnet:     a struct with the following options to set:
%       .nickname       = nickname of the sub-analysis
%       .netmetric      = 'nbs' or 'maxT'
%       .fb_interest    = 'alpha', 'beta', 'theta', 'Lgamma', 'Hgamma'
%       .edge_weight    = 1     if you wanted edge thickness by diff in
%                               connectivity value between the two groups;
%                         0     for all edges having the same thickness
%       .colour         = 1     for one colour for all nodes
%                         2     for colour range by number of connections 
%                               per node (node degree)
%       .size           = 1     for size of each node is the same
%                         2 	for size based on group connectivity
%                         3     for size based on group connectivity difference
%       .grp            = division of participants by index into groups,
%                         with a different cell entry per group (e.g.,
%                         brainnet.grp = {1:19, 20:32}.
%       .noi            = nodes of interest; if left empty, all nodes with
%                         some connection to another node will be included 
%                         in visualization.
% 
% OUTPUT-------------------------------------------------------------------
% [nickname]_[fb_interest].edge:   design matrix; if you have a number of
%                                  nodes of interest (NOI), then this is a (NOI x NOI) text
%                                  file where each entry represents the edge weight between
%                                  the node row-column combination. Only include the nodes you
%                                  want to visualize. 
% [nickname]_[fb_interest].node:   info on nodes in a (NOI x 6 col) space-delimited text file
%                                  for each row (i.e., each node of interest):
%                                   - col 1-3 = coordinates of the node
%                                   - col 4 = sphere size
%                                   - col 5 = sphere colour
%                                   - col 6 = region label


% IMPORTANT: 
%   CURRENTLY MAINLY WORKS FOR 'NBS'
%   NOTE FOR PLS:
%      - positive (_pos.node/.edge) files refer to POSITIVE saliences not 
%        necessarily an INCREASED connectivity network (looking at your 
%        contrast from PLS)
%      - negative (_neg.node/.edge) files refer to NEGATIVE saliences not 
%        necessarily a DECREASED connectivity network (looking at your 
%        contrast from PLS)
% 
% Written by Sonya, edited by Liz, function-ified by Julie (March 2020)
    
%% setup

if strcmp(brainnet.netmetric, 'nbs')

    % load nbs results file: */analysis/group/capybara_NBSResults_alpha.mat
    load([paths.anout_grp '/' brainnet.nickname '_NBSResults_' brainnet.fb_interest '.mat'])

    % load the NBS data matrix: */analysis/group/capybara_datamatrix_alpha.mat
    load([paths.anout_grp '/' brainnet.nickname '_datamatrix_' brainnet.fb_interest '.mat'])
    
elseif strcmp(brainnet.netmetric, 'maxT')
    load([paths.anout '/NBS_treatmenttype/NBS' brainnet.fb_interest '_Hz.mat']) % load data corresponding to freq band
end

% get difference between groups
diff_groupadjmat = mean(data_matrix(:,:,brainnet.grp{1}),3) - ...
                    mean(data_matrix(:,:,brainnet.grp{2}),3);

%% EDGE FILE

% prepare binary mask (sig_nodes)
if strcmp(brainnet.netmetric,'nbs')
    sig_nodes   = full(nbs.NBS.con_mat{1,1}+nbs.NBS.con_mat{1,1}');
elseif strcmp(brainnet.netmetric,'maxT')
    sig_nodes   = zeros(90,90);
    for i = 1:length(brainnet.maxT_NOI)
        if length(brainnet.maxT_NOI{i}) > 1
            sig_nodes(brainnet.maxT_NOI{i}(1), brainnet.maxT_NOI{i}(2:end)) = 1;
        end
    end
    sig_nodes = sig_nodes + sig_nodes';
end

% isolate only nodes we want to look at
if isempty(brainnet.noi) % if none specified, then use all nodes with existing connections
    noi = sum(sig_nodes) > 0;
else
    noi = brainnet.noi;
end

if brainnet.edge_weight == 1
    sig_nodes = sig_nodes .* diff_groupadjmat;
end

if strcmp(brainnet.netmetric,'nbs')
    dlmwrite([paths.anout_grp '/' brainnet.nickname '_' brainnet.fb_interest '.edge'],sig_nodes(noi,noi),'delimiter',' ','precision','%d');
elseif strcmp(brainnet.netmetric,'maxT')
    dlmwrite([paths.anout '/NBS_treatmenttype/NBS' brainnet.fb_interest '_maxT.edge'],sig_nodes(noi,noi),'delimiter',' ','precision','%d');
end

%% NODE FILE

% SET COLOUR---------------------------------------------------------------
if brainnet.colour == 1
    node_colour                     = squeeze(sum(sig_nodes,1));
    node_colour(node_colour > 0)    = 1;
elseif brainnet.colour == 2
    node_colour                     = squeeze(nansum(sig_nodes(noi,noi),1));
elseif brainnet.colour==3
    error('Sorry - brainnet.colour = 3 is not yet a working option!')
%     % define lobes
%     lobes                           = zeros(1,90);
%     lobes([1:28,69:70])             = 1; % F = frontal
%     lobes([29:36])                  = 2; % I = insular and cingulate gyri
%     lobes([37:42,55:56,79:90])      = 3; % T = temporal
%     lobes([43:54])                  = 4; % O = occipital
%     lobes([57:68])                  = 5; % P = parietal
%     lobes([71:78])                  = 6; % S = subcortical
%     % lobes([91:116])               = 7; % C = cerebellar regions
%     if strcmp(brainnet.netmetric,'nbs')
%         lobes                       = repmat(lobes',[1,2]);
%         define colours
%         node_colour                 = squeeze(sum(sig_nodes,1));
%         node_colour(node_colour > 0)= 1;
%         node_colour                 = node_colour.*lobes(1:size(sig_nodes,1),:);
%         node_colour                 = lobes(1:size(sig_nodes,1),:);
%     elseif strcmp(brainnet.netmetric,'nbs')
%         lobes                       = lobes';
%         % define colours
%         node_colour                 = squeeze(sum(sig_nodes,1));
%         node_colour(node_colour > 0)= 1;
%         node_colour                 = node_colour.*lobes(1:size(sig_nodes,1));
%         node_colour                 = lobes(1:size(sig_nodes,1)); 
%     end
end

% NODE SIZE----------------------------------------------------------------
if strcmp(brainnet.netmetric,'nbs') || strcmp(brainnet.netmetric, 'maxT')
    if brainnet.size_nodes == 1 % size of each node identical
        node_size                   = ones(1,size(sig_nodes(noi,noi),1));
    elseif brainnet.size_nodes == 2 % size based on group connectivity
        group_nodes                 = groupadjmat(:,:,brainnet.fb_interest,1).*sig_nodes; 
        node_size                   = strengths_und(group_nodes);
    elseif brainnet.size_nodes == 3 % size based on difference between groups
        diff_nodes                  = diff_groupadjmat.*sig_nodes(noi,noi);
        diff_nodes(logical(eye(sum(noi)))) = 0;
        node_size                   = strengths_und(abs(diff_nodes));
    end
end
node_size=abs(round(node_size',4)); % gets rid of negatives... so size is accurately depicted

% ASSEMBLE INTO A NODE FILE------------------------------------------------
if strcmp(brainnet.netmetric,'nbs')
    region_labels                   = nbs.NBS.node_label(noi);
    node_file                       = cat(2,nbs.NBS.node_coor(noi,:),node_colour');
    node_file                       = cat(2,node_file,node_size); 
    node_file                       = cat(2,num2cell(node_file),region_labels);
    % save .node file (convert from cell array to string array in order to save
    fid                             = fopen([paths.anout_grp '/' brainnet.nickname '_' brainnet.fb_interest '.node'],'wt');
    for j = 1:size(node_file,1)
        tmp_name                    = strsplit(region_labels{j,1},'_');
        new_region_label            = strjoin(tmp_name,'.');
        fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\n',node_file{j,1},node_file{j,2},node_file{j,3},node_file{j,4},node_file{j,5},new_region_label);
    end
    fclose(fid);
elseif strcmp(brainnet.netmetric, 'maxT')
    % load sample nbs
    load(['/mnt/sda/juanita/MEGneto/analysis/left/analysis/NBS_RAD_vs_TDC/RAD_vs_TDC_Lgamma_extent.mat'])
    
    region_labels                   = nbs.NBS.node_label(noi);
    node_file                       = cat(2,nbs.NBS.node_coor(noi,:),node_colour');
    node_file                       = cat(2,node_file,node_size); 
    node_file                       = cat(2,num2cell(node_file),region_labels);
    % save .node file (convert from cell array to string array in order to save)
    fid                             = fopen([paths.anout '/NBS_treatmenttype/NBS' brainnet.fb_interest '_maxT.node'],'wt');
    for j = 1:size(node_file,1)
        tmp_name                    = strsplit(region_labels{j,1},'_');
        new_region_label            = strjoin(tmp_name,'.');
        fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\n',node_file{j,1},node_file{j,2},node_file{j,3},node_file{j,4},node_file{j,5},new_region_label);
    end
    fclose(fid);
end

%% Say goodbye!
disp('Done making BNV *.edge and *.node files!');