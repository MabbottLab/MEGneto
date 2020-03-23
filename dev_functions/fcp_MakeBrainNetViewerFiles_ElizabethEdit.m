% Makes .node and .edge files for BrainNet Viewer
%       NOTE FOR PLS:
%               -positive (_pos.node/.edge) files refer to POSITIVE saliences not necessarily an INCREASED connectivity network
%                       (looking at your contrast from PLS)
%               -negative (_neg.node/.edge) files refer to NEGATIVE saliences not necessarily a DECREASED connectivity network
%                       (looking at your contrast from PLS)
% clear;
% clc;
% informationENT VERSION: 
% brainnet.pls_path = '/data/diana/exercise/ANALYSIS/diff_ICA/PLS_2-3Hz_lv1.mat';
% brainnet.pls_path = '/home/sonya/Data_exercise/results_dAttention_nocrossover6/PLS_60-100Hz_lv1.mat';
clear node_size;
addpath('/home/liz/matlab/BCT/2016_01_16_BCT');
addpath('/home/liz/matlab/BrainNetViewer')
brainnet.nbs_path = ' ';
brainnet.netmetric = 'nbs'; % network metric ('pls' or 'nbs');
brainnet.connmetric = 'pli'; % connectivity metric
brainnet.fb_interest= 4; % frequency band of interest (integer refering to frequency band present in p.filt_freqs) 

% output directory
brainnet.ouput_dir = '/mnt/sda/juanita';

% output file name (excluding .node and .edge)
brainnet.file_name= 'RADvsSURG';

% colour of nodes (determines colour of nodes in BrainNet Viewer)
%[1 = one colour for all nodes]
%[2 = colour range by number of connections per node (node degree)]
%[3 = colour range by lobes]
brainnet.colour = 2;

% size of nodes (determines size of nodes in BrainNet Viewer)
%[1 = size of each node is the same]
%[2 = size based on group connectivity]
%[3 = size based on group connectivity difference]
%[4 = size based on saliences - how much each region contributes to the LV (PLS only)]
brainnet.size_nodes = 3;

% p-structures to load (need if size_nodes = 2 or 3)
% brainnet.group_paths{1} = '/data/diana/exercise/pre/groupANALYSIS/pINFO.mat'; % (need if size_nodes = 2 or 3)
% brainnet.group_paths{2} = '/data/diana/exercise/post/groupANALYSIS/pINFO.mat'; % (need if size_nodes = 3)
% brainnet.group_paths{1} = '/home/sonya/Data_exercise/T1/groupANALYSIS/pINFO.mat';
% brainnet.group_paths{2} = '/home/sonya/Data_exercise/T2/groupANALYSIS/pINFO.mat';


%% 

% load needed files
% file_in='cerebellum_region_coords.mat';
file_in='/mnt/sda/juanita/MEGneto_old/functions/functions/AAL_Atlas/AAL116_region_coords.mat'; % AAL_90sources_region_coords
load(file_in);
% file_in='cerebellum_region_names.mat';
file_in='/mnt/sda/juanita/MEGneto_old/functions/functions/AAL_Atlas/AAL116_region_names.mat';% AAL_90sources_region_labels
load(file_in);

% if strcmp(brainnet.netmetric,'pls') 
%     load(brainnet.pls_path);
% elseif strcmp(brainnet.netmetric,'nbs') 
%     load(brainnet.nbs_path);
% end

% if brainnet.size_nodes == 2 ||  brainnet.size_nodes == 3
%     % load p-structures
%     numGroups = length(brainnet.group_paths);
%     for i=1:numGroups
%         % load p-structure
%         load(brainnet.group_paths{i},'p');
%         % %%% INPUTS %%% (update anonymous function's workspace if any edits to subject paths made)
%         p.paths.connmat =           @(id,postfix) [p.paths.output_dir(id), p.subj.ID{id}, '_', postfix,'.mat'];
%         brainnet.p{i} = p;
%     end
% 
%     % load data
%     matlinks = {};
%     for gg = 1:numGroups
%         numSubj = length(brainnet.p{gg}.subj.ID);
%         for ss = 1:numSubj
%             matlinks{gg}{ss} = matfile(brainnet.p{gg}.paths.connmat(ss,brainnet.connmetric));                
%         end
%     end
% 
%     % get dimensions and frequency bands
%     m = matfile(brainnet.p{gg}.paths.connmat(1,brainnet.connmetric));
%     [numSources, ~, ~, numFreqs] = size(m.adjmat); % must be the same for each subject
%     filt_freqs=brainnet.p{1}.filt_freqs;
% 
%     % Make group averages and differences
%     for gg = 1:numGroups
%         numSubj = length(brainnet.p{gg}.subj.ID); % can differ per group
%         fprintf('Making averages...\n');
%         for ss = 1:numSubj
%             temp_subjadjmat=matlinks{1, gg}{1, ss}.adjmat;
%             [~, ~, numTrials, ~] = size(temp_subjadjmat); % can differ per subject
%             % set diagonal to NaN, so that the mean isn't biased.
%             eyenan = ones(numSources);
%             eyenan(eye(numSources)==1) = NaN;
%             temp_subjadjmat=temp_subjadjmat.*repmat(eyenan,1,1,numTrials,numFreqs);
%             % Average Trials [nchan nchan ntrial nfreq] -> [nchan nchan nfreq] 
%             temp_subjadjmat=squeeze(mean(temp_subjadjmat,3,'omitnan'));
%             groupadjmat(:,:,:,ss,gg)=temp_subjadjmat;  % [nchan nchan nfreq nsubj ngroup]
%         end
%     end
%     % Rearrange matrix [nchan nchan nsubj nfreq ngroup]
%     groupadjmat=permute(groupadjmat,[1,2,4,3,5]); 
%     % Average Subj [nchan nchan nsubj nfreq ngroup] -> [nchan nchan nfreq ngroup]
%     groupadjmat=squeeze(mean(groupadjmat,3,'omitnan'));
%     % set diagonal to 1 (fully correlated)
%     groupadjmat(logical(repmat(eye(numSources),1,1,numFreqs,numGroups))) = 1;
% 
%     brainnet = rmfield(brainnet,{'p'});
%     clear p matlinks m 
% 
%     % Group differences
%     numDiff= numGroups-1; % should be = 1
%     if numDiff~=0
%         diff_groupadjmat=zeros(numSources,numSources,numFreqs,numDiff);
%         for dd = 1:numDiff
%             diff_groupadjmat(:,:,:,dd)=groupadjmat(:,:,:,dd+1)-groupadjmat(:,:,:,dd);
%         end
%     end
% end

%%
% THRESHOLD SALIENCES
thresh_percent = 0.52;
%%%%%%%%%%%%%
temp_thresholded_saliences = nan(size(thresholded_saliences));
temp_thresholded_saliences(thresholded_saliences > 0) = thresholded_saliences(thresholded_saliences > 0);
salience_min = min(temp_thresholded_saliences(:),[],'omitnan');
salience_max = max(temp_thresholded_saliences(:),[],'omitnan');
pos_salience_threshold = ((salience_max-salience_min)*thresh_percent) + salience_min;

temp_thresholded_saliences = nan(size(thresholded_saliences));
temp_thresholded_saliences(thresholded_saliences < 0) = thresholded_saliences(thresholded_saliences < 0);
salience_max = min(temp_thresholded_saliences(:),[],'omitnan');
salience_min = max(temp_thresholded_saliences(:),[],'omitnan');
neg_salience_threshold = salience_max - ((salience_max-salience_min)*thresh_percent);

thresholded_saliences(thresholded_saliences>neg_salience_threshold & thresholded_saliences<pos_salience_threshold) = 0;
%%%%%%%%%%%%%%

%%
% prepare binary mask (sig_nodes)
if strcmp(brainnet.netmetric,'pls')
    sig_nodes = cat(3,thresholded_saliences > 0,thresholded_saliences < 0);
elseif strcmp(brainnet.netmetric,'nbs')
    sig_nodes=full(nbs.NBS.con_mat{1,1}+nbs.NBS.con_mat{1,1}');
end

% output file
file_out=strcat(brainnet.ouput_dir,'/',brainnet.file_name);
% 
% save .edge file
if strcmp(brainnet.netmetric,'pls')
    dlmwrite([file_out,'_pos.edge'],sig_nodes(:,:,1),'delimiter',' ','precision','%d');
    dlmwrite([file_out,'_neg.edge'],sig_nodes(:,:,2),'delimiter',' ','precision','%d');
elseif strcmp(brainnet.netmetric,'nbs')
    dlmwrite([file_out,'.edge'],sig_nodes,'delimiter',' ','precision','%d');
end

%%
% thresholded edge file: make edge file for specific AAL region indicated
% used to visualize connections with hub nodes only 

% % OLD:
sig_nodes_temp = sig_nodes;
sig_nodes_temp(:,[1:18],:) =0;
% 
%% to visualize all connections with 1 node of interest
% sig_nodes_temp = sig_nodes;
% sig_nodes_temp([1,3:116],[1,3:116],:) =0; % nodes NOT visualized -> indicated by 0

dlmwrite([file_out,'_RADSURG_pos.edge'],sig_nodes_temp(:,:,1),'delimiter',' ','precision','%d');
dlmwrite([file_out,'_RADSURG_neg.edge'],sig_nodes_temp(:,:,2),'delimiter',' ','precision','%d');

% %to visualize connections with all nodes in lobe of interest 
% sig_nodes_temp = sig_nodes;
% sig_nodes_temp([91:116],[91:116])=0; % use brainnet colour lobes to indicate nodes NOT to visualize 
% %%(indicate nodes not in lobe of interest)
% %%example: want to visualize FRONTAL LOBE = 1:28,69:70 (nodes 1-28,69,70) - therefore indicate set 29-68, 71-116=0 (edges=no, not displayed)
% dlmwrite([file_out,'_cerebellum_pos.edge'],sig_nodes_temp(:,:,1),'delimiter',' ','precision','%d');
% dlmwrite([file_out,'_cerebellum_neg.edge'],sig_nodes_temp(:,:,2),'delimiter',' ','precision','%d');

%%
% node colour
if brainnet.colour==1
    node_colour=squeeze(sum(sig_nodes,1));
    node_colour(node_colour > 0) = 1;
elseif brainnet.colour==2
    node_colour=squeeze(sum(sig_nodes,1));
elseif brainnet.colour==3
    % define lobes
    lobes=zeros(1,116);
    lobes([1:28,69:70])=1; % F = frontal
    lobes([29:36])=2; % I = insular and cingulate gyri
    lobes([37:42,55:56,79:90])=3; % T = temporal
    lobes([43:54])=4; % O = occipital
    lobes([57:68])=5; % P = parietal
    lobes([71:78])=6;% S = subcortical
    lobes([91:116])=7; % C = cerebellar regions
    if strcmp(brainnet.netmetric,'nbs')
        lobes=repmat(lobes',[1,2]);
        % define colours
%         node_colour=squeeze(sum(sig_nodes,1));
%         node_colour(node_colour > 0) = 1;
%         node_colour=node_colour.*lobes(1:size(sig_nodes,1),:);
        node_colour = lobes(1:size(sig_nodes,1),:);
%     elseif strcmp(brainnet.netmetric,'nbs')
%         lobes=lobes';
%         % define colours
% %         node_colour=squeeze(sum(sig_nodes,1));
% %         node_colour(node_colour > 0) = 1;
% %         node_colour=node_colour.*lobes(1:size(sig_nodes,1));
%         node_colour=lobes(1:size(sig_nodes,1)); 
    end
end
%to only visualize edges for specific node colours
% for lobes = 1:1
%     sig_nodes_temp = sig_nodes;
%     sig_nodes([1:28,69:70],[1:28,69:70],:) = 0;
% end
%    
% 
% dlmwrite([file_out,'_posTEST.edge'],sig_nodes_temp(:,:,1),'delimiter',' ','precision','%d');
% dlmwrite([file_out,'_negTEST.edge'],sig_nodes_temp(:,:,2),'delimiter',' ','precision','%d');

% node size
if strcmp(brainnet.netmetric,'nbs')
    if brainnet.size_nodes==1 % size of each node identical
        node_size= ones(2,size(sig_nodes,1));
    elseif brainnet.size_nodes==2 % size based on group connectivity
        group_nodes(:,:,1)=groupadjmat(:,:,brainnet.fb_interest,1).*sig_nodes(:,:,1); % pos binary mask for nodes in network
        group_nodes(:,:,2)=groupadjmat(:,:,brainnet.fb_interest,1).*sig_nodes(:,:,2); % neg binary mask for nodes in network
        node_size(1,:) = strengths_und(group_nodes(:,:,1));
        node_size(2,:) = strengths_und(group_nodes(:,:,2));
    elseif brainnet.size_nodes==3 % size based on difference between groups
        diff_nodes(:,:,1)=diff_groupadjmat(:,:,brainnet.fb_interest,1).*sig_nodes(:,:,1); % pos binary mask for nodes in network
        diff_nodes(:,:,2)=diff_groupadjmat(:,:,brainnet.fb_interest,1).*sig_nodes(:,:,2); % neg binary mask for nodes in network
        node_size(1,:) = strengths_und(diff_nodes(:,:,1));
        node_size(2,:) = strengths_und(diff_nodes(:,:,2));
    elseif brainnet.size_nodes==4 % size based on saliences
%         node_size(1,:) = strengths_und(thresholded_saliences.*(thresholded_saliences > 0)); % pos binary mask for nodes in network
%         node_size(2,:) = strengths_und(thresholded_saliences.*(thresholded_saliences < 0)); % neg binary mask for nodes in network
        %node_size 116x2
        mask_pos = pos_salience_strength;
        mask_neg = neg_salience_strength;
%         mask_pos(mask_pos < 0.1)=0;
%         mask_neg(mask_neg > -0.1)=0;
        node_size(1,:) = mask_pos; 
        node_size(2,:) = mask_neg;
    end
% elseif strcmp(brainnet.netmetric,'nbs')
%     if brainnet.size_nodes==1 % size of each node identical
%         node_size= ones(1,size(sig_nodes,1));
%     elseif brainnet.size_nodes==2 % size based on group connectivity
%         group_nodes=groupadjmat(:,:,brainnet.fb_interest,1).*sig_nodes; 
%         node_size = strengths_und(group_nodes);
%     elseif brainnet.size_nodes==3 % size based on difference between groups
%         diff_nodes=diff_groupadjmat(:,:,brainnet.fb_interest,1).*sig_nodes;
%         node_size = strengths_und(diff_nodes);
%     end
end
node_size=abs(round(node_size',4)); % gets rid of negatives... so size is accurately depicted

% make .node file
if strcmp(brainnet.netmetric,'pls')
    for i=1:2
        node_file=cat(2,region_coords,node_colour(:,i));
        node_file=cat(2,node_file,node_size(:,i)); 
        node_file=cat(2,num2cell(node_file),region_labels);
        % save .node file (convert from cell array to string array in order to save)
        if i==1
            fid=fopen([file_out,'_pos.node'],'wt');
        else
            fid=fopen([file_out,'_neg.node'],'wt');
        end 
        for j=1:size(node_file,1);
            tmp_name = strsplit(region_labels{j,1},'_');
            new_region_label= strjoin(tmp_name,'.');
            fprintf(fid,'%d\t%d\t%d\t%d\t%f\t%s\n',node_file{j,1},node_file{j,2},node_file{j,3},node_file{j,4},node_file{j,5},new_region_label);
        end
        fclose(fid);
    end
% elseif strcmp(brainnet.netmetric,'nbs')
%     node_file=cat(2,region_coords,node_colour);
%     node_file=cat(2,node_file,node_size); 
%     node_file=cat(2,num2cell(node_file),region_labels);
%     % save .node file (convert from cell array to string array in order to save)
%     fid=fopen([file_out,'.node'],'wt');
%     for j=1:size(node_file,1);
%         tmp_name = strsplit(region_labels{j,1},'_');
%         new_region_label= strjoin(tmp_name,'.');
%         fprintf(fid,'%d\t%d\t%d\t%d\t%f\t%s\n',node_file{j,1},node_file{j,2},node_file{j,3},node_file{j,4},node_file{j,5},new_region_label);
%     end
%     fclose(fid);
end

disp('Done');