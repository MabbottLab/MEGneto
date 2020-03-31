function make_BNV_ready(paths, brainnet)

% To make *.node and *.edge files for viewing connectivity results from PLS
% or NBS on BrainNet Viewer (BNV). 
% 
% INPUT--------------------------------------------------------------------
% paths:        the paths struct, as we always input
% brainnet:     a struct with the following options to set:
%       .netmetric      = 'pls' or 'nbs'
%       .connmetric     = 'pli' or 'wpli' or 'wpli_debiased' ... etc.
%       .fb_interest    = 'alpha', 'beta', 'theta', 'Lgamma', 'Hgamma'
%       .file_name      = make this the same as the folder holding the NBS
%                         results and data matrices from make_NBS_ready.m
%       .edge_weight    = 1     if you wanted edge thickness by diff in
%                               connectivity value between the two groups;
%                         0     for all edges having the same thickness
%       .colour         = 1     for one colour for all nodes
%                         2     for colour range by number of connections 
%                               per node (node degree)
%                         3     for colour range by lobes
%       .size           = 1     for size of each node is the same
%                         2 	for size based on group connectivity
%                         3     for size based on group connectivity difference
%                         4     for size based on saliences - how much each 
%                               region contributes to the LV (PLS only)]
%       .grp            = division of participants by index into groups,
%                         with a different cell entry per group (e.g.,
%                         brainnet.grp = {1:19, 20:32}.
%       .noi            = nodes of interest; if left empty, all nodes with
%                         some connection to another node will be included 
%                         in visualization.
% 
% OUTPUT-------------------------------------------------------------------
% [fb_interest].edge:   design matrix; if you have a number of
%                       nodes of interest (NOI), then this is a (NOI x NOI) text
%                       file where each entry represents the edge weight between
%                       the node row-column combination. Only include the nodes you
%                       want to visualize. 
% [fb_interest].node:   info on nodes in a (NOI x 6 col) space-delimited text file
%                       for each row (i.e., each node of interest):
%                           - col 1-3 = coordinates of the node
%                           - col 4 = sphere size
%                           - col 5 = sphere colour
%                           - col 6 = region label

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

% output directory
brainnet.ouput_dir      = [paths.anout '/' brainnet.file_name '/'];

% load nbs results file
load(char(glob([paths.anout '/' brainnet.file_name '/*' brainnet.fb_interest '_' brainnet.nbs_type '*'])))

% load the NBS data matrix
load(char(glob([paths.anout '/' brainnet.file_name '/*' brainnet.fb_interest '_Hz.mat'])))

% get difference between groups
diff_groupadjmat = mean(data_matrix(:,:,brainnet.grp{1}),3) - ...
                    mean(data_matrix(:,:,brainnet.grp{2}),3);

%% EDGE FILE

% prepare binary mask (sig_nodes)
if strcmp(brainnet.netmetric,'pls')
    sig_nodes   = cat(3,thresholded_saliences > 0,thresholded_saliences < 0);
elseif strcmp(brainnet.netmetric,'nbs')
    sig_nodes   = full(nbs.NBS.con_mat{1,1}+nbs.NBS.con_mat{1,1}');
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

if strcmp(brainnet.netmetric,'pls')
    dlmwrite([brainnet.ouput_dir brainnet.fb_interest '_pos.edge'],sig_nodes(:,:,1),'delimiter',' ','precision','%d');
    dlmwrite([brainnet.ouput_dir brainnet.fb_interest '_neg.edge'],sig_nodes(:,:,2),'delimiter',' ','precision','%d');
elseif strcmp(brainnet.netmetric,'nbs')
    dlmwrite([brainnet.ouput_dir brainnet.fb_interest '.edge'],sig_nodes(noi,noi),'delimiter',' ','precision','%d');
end

%% NODE FILE

% SET COLOUR---------------------------------------------------------------
if brainnet.colour == 1
    node_colour                     = squeeze(sum(sig_nodes,1));
    node_colour(node_colour > 0)    = 1;
elseif brainnet.colour == 2
    node_colour                     = squeeze(nansum(sig_nodes(noi,noi),1));
elseif brainnet.colour==3
    % define lobes
    lobes                           = zeros(1,90);
    lobes([1:28,69:70])             = 1; % F = frontal
    lobes([29:36])                  = 2; % I = insular and cingulate gyri
    lobes([37:42,55:56,79:90])      = 3; % T = temporal
    lobes([43:54])                  = 4; % O = occipital
    lobes([57:68])                  = 5; % P = parietal
    lobes([71:78])                  = 6; % S = subcortical
    % lobes([91:116])               = 7; % C = cerebellar regions
    if strcmp(brainnet.netmetric,'nbs')
        lobes                       = repmat(lobes',[1,2]);
        define colours
        node_colour                 = squeeze(sum(sig_nodes,1));
        node_colour(node_colour > 0)= 1;
        node_colour                 = node_colour.*lobes(1:size(sig_nodes,1),:);
        node_colour                 = lobes(1:size(sig_nodes,1),:);
    elseif strcmp(brainnet.netmetric,'nbs')
        lobes                       = lobes';
        % define colours
        node_colour                 = squeeze(sum(sig_nodes,1));
        node_colour(node_colour > 0)= 1;
        node_colour                 = node_colour.*lobes(1:size(sig_nodes,1));
        node_colour                 = lobes(1:size(sig_nodes,1)); 
    end
end

% NODE SIZE----------------------------------------------------------------
if strcmp(brainnet.netmetric,'pls')
    if brainnet.size_nodes == 1 % size of each node identical
        node_size                   = ones(2,size(sig_nodes(noi,noi),1));
    elseif brainnet.size_nodes == 2 % size based on group connectivity
        group_nodes(:,:,1)          = groupadjmat(:,:,brainnet.fb_interest,1).*sig_nodes(:,:,1); % pos binary mask for nodes in network
        group_nodes(:,:,2)          = groupadjmat(:,:,brainnet.fb_interest,1).*sig_nodes(:,:,2); % neg binary mask for nodes in network
        node_size(1,:)              = strengths_und(group_nodes(:,:,1));
        node_size(2,:)              = strengths_und(group_nodes(:,:,2));
    elseif brainnet.size_nodes == 3 % size based on difference between groups
        diff_nodes(:,:,1)           = diff_groupadjmat(:,:,brainnet.fb_interest,1).*sig_nodes(:,:,1); % pos binary mask for nodes in network
        diff_nodes(:,:,2)           = diff_groupadjmat(:,:,brainnet.fb_interest,1).*sig_nodes(:,:,2); % neg binary mask for nodes in network
        node_size(1,:)              = strengths_und(diff_nodes(:,:,1));
        node_size(2,:)              = strengths_und(diff_nodes(:,:,2));
    elseif brainnet.size_nodes == 4 % size based on saliences
%       node_size(1,:)              = strengths_und(thresholded_saliences.*(thresholded_saliences > 0)); % pos binary mask for nodes in network
%       node_size(2,:)              = strengths_und(thresholded_saliences.*(thresholded_saliences < 0)); % neg binary mask for nodes in network
%       node_size 116x2
        mask_pos                    = pos_salience_strength;
        mask_neg                    = neg_salience_strength;
%       mask_pos(mask_pos < 0.1)    = 0;
%       mask_neg(mask_neg > -0.1)   = 0;
        node_size(1,:)              = mask_pos; 
        node_size(2,:)              = mask_neg;
    end
elseif strcmp(brainnet.netmetric,'nbs')
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
if strcmp(brainnet.netmetric,'pls')
    for i = 1:2
        node_file                   = cat(2,region_coords,node_colour(:,i));
        node_file                   = cat(2,node_file,node_size(:,i)); 
        node_file                   = cat(2,num2cell(node_file),region_labels);
        % save .node file (convert from cell array to string array in order to save)
        if i == 1
            fid                     = fopen([file_out,'_pos.node'],'wt');
        else
            fid                     = fopen([file_out,'_neg.node'],'wt');
        end 
        for j = 1:size(node_file,1)
            tmp_name                = strsplit(region_labels{j,1},'_');
            new_region_label        = strjoin(tmp_name,'.');
            fprintf(fid,'%d\t%d\t%d\t%d\t%f\t%s\n',node_file{j,1},node_file{j,2},node_file{j,3},node_file{j,4},node_file{j,5},new_region_label);
        end
        fclose(fid);
    end
elseif strcmp(brainnet.netmetric,'nbs')
    region_labels                   = nbs.NBS.node_label(noi);
    node_file                       = cat(2,nbs.NBS.node_coor(noi,:),node_colour');
    node_file                       = cat(2,node_file,node_size); 
    node_file                       = cat(2,num2cell(node_file),region_labels);
    % save .node file (convert from cell array to string array in order to save)
    fid                             = fopen([brainnet.ouput_dir brainnet.fb_interest '.node'],'wt');
    for j = 1:size(node_file,1)
        tmp_name                    = strsplit(region_labels{j,1},'_');
        new_region_label            = strjoin(tmp_name,'.');
        fprintf(fid,'%d\t%d\t%d\t%d\t%f\t%s\n',node_file{j,1},node_file{j,2},node_file{j,3},node_file{j,4},node_file{j,5},new_region_label);
    end
    fclose(fid);
end

%% Say goodbye!
disp('Done making BNV *.edge and *.node files!');
