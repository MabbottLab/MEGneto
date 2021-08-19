%% SCRAP CODE
% This script contains several scrap bits of code that may come in handy
% and are left here for future convenience. This code was not used for
% analysis.

%% See interpolation on dipole grid

% load template grid with 5mm resolution
% template_grid       = load([fieldtrip_path '/template/sourcemodel/standard_sourcemodel3d5mm.mat']);
% template_grid       = template_grid.sourcemodel;

% interpolate atlas to template grid
% cfg                 = [];
% cfg.interpmethod    = 'nearest';
% cfg.parameter       = 'tissue';
% sourcemodel_mmp     = ft_sourceinterpolate(cfg, atlas, template_grid);

% for running a big loop
of_interest = [1, 10, 11, 18]; % 1 = V1, 10 = FEF, 11 = PEF, 18 = FFC

% for regions found in the atlas
idx = 18;
roi = find(source_atlas.tissue == idx);
roi_pos = source_atlas.pos(roi, :);

% for external coordinates
roi_pos = [-8 -76 10]; % .e.g, L-V1 MNI coordinates

ft_plot_mesh(source_atlas.pos(template_grid.inside,:), 'vertexsize', 3); % source_atlas = template_grid
hold on;
for i = 1:length(roi)
    ft_plot_mesh(source_atlas.pos(roi(i), :), 'vertexcolor', 'r', 'vertexsize', 20); 
    disp(source_atlas.pos(roi(i), :))
end
title(source_atlas.tissuelabel{idx}, 'interpreter', 'none');
hold off 
view([-270 0]);




var_explained  = NaN(1, ... % set up empty matrix to store variance explained of first principial component for each trial and region of interest 
                         length(projection.trial), ...
                         length(source_atlas.tissuelabel));

%% Plotting 12 V1 ROIs for 1 trial and 1 ppt
catmatrix = NaN(12,...
                length(projection.time), ... % set up empty matrix to store reconstructed timeseries for each trial and region of interest
                length(projection.trial), ...
                length(source_atlas.tissuelabel));   % over all ROI timeseries across trial                     
for t = 1:projection.df
    %%% AND FOR EACH NODE ---------------------------------------------
    for i = 1
        % identify source coords that fall within ROI
        node                     = find(source_atlas.tissue == i);
        source_timeseries        = cell2mat(projection.trial(t).mom(node)); % get the timeseries; num_nodes x time
        % ori_region               = cell2mat(projection.trial(t).ori(node)); % orientations; num_nodes x time
        for j = 1:size(source_timeseries, 1)
            catmatrix(j,:,t,i) = source_timeseries(j,:);
        end
    end
end

for i = 1:12
    subplot(12,1,i)
    plot(projection.time(1:1000), catmatrix(i,1:1000,1,1))
end

%% Plotting V1 average across 12 regions, across all trials
 for t = 1:projection.df
     %%% AND FOR EACH NODE ---------------------------------------------
     for i = 18
         % identify source coords that fall within ROI
         node                     = find(source_atlas.tissue==i);
         source_timeseries        = cell2mat(projection.trial(t).mom(node)); % get the timeseries; num_nodes x time
         % ori_region               = cell2mat(projection.trial(t).ori(node)); % orientations; num_nodes x time
         
         % IF NODE EXISTS
         if size(source_timeseries, 1) >= 1
             if config.beamforming.rep_timeseries == "mean"
                 catmatrix(:,t,i) = nanmean(source_timeseries,1); % take avg across source points
             elseif config.beamforming.rep_timeseries == "pca"
                 [~, score, ~, ~, explained] = pca(transpose(source_timeseries)); % perform pca
                 catmatrix(:,t,i) = transpose(score(:, 1)); % store first principal component across timeseries
                 var_explained(:,t,i) = explained(1);
             end
         else
             warning('No sources in ROI %s.\n',atlas.tissuelabel{i});
         end
     end
 end
trials_avg = nanmean(catmatrix, 2);
plot(projection.time(1:1000), trials_avg(1:1000,1,i))
title('FFA (18 in Glasser): average across dipoles and average across trials')

%% Plotting individuals nodes of V1 (12 total) but averaged across 82 trials
catmatrix = NaN(12,...
                length(projection.time), ... % set up empty matrix to store reconstructed timeseries for each trial and region of interest
                length(projection.trial), ...
                length(source_atlas.tissuelabel));   % over all ROI timeseries across trial                     
for t = 1:projection.df
    %%% AND FOR EACH NODE ---------------------------------------------
    for i = 1
        % identify source coords that fall within ROI
        node                     = find(source_atlas.tissue == i);
        source_timeseries        = cell2mat(projection.trial(t).mom(node)); % get the timeseries; num_nodes x time
        % ori_region               = cell2mat(projection.trial(t).ori(node)); % orientations; num_nodes x time
        for j = 1:size(source_timeseries, 1)
            catmatrix(j,:,t,i) = source_timeseries(j,:);
        end
    end
end

trials_avg_per_node = nanmean(catmatrix, 3);
for i = 1:12
    subplot(12,1,i)
    plot(projection.time(1:1000), trials_avg_per_node(i,1:1000,1,1))
end

%% START OF AVG ACROSS PPTS AND TRIALS 
%% Plotting across participants and across trials with Glasser atlas
config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp4';

% check for matched MRI and MEG data
subj_match  = freeviewing_ds_pid_match(paths,step);
ssSubjPath  = @(x) paths.(subj_match.pid{x});

rangeOFsubj = 1:length(subj_match.pid);

all_catmatrix = NaN(2401,...
                    130,...
                    360,...
                    length(subj_match.pid));

for ss = rangeOFsubj % for each participant that has matched MEG/MRI data
    prior_dir = pwd();
    cd(['/home/dmatic/MEGProjects/analysis/OIRM_test_May25/analysis/' subj_match.pid{ss}])
    current_catmatrix = load('atlas_beamforming_results.mat');
    current_catmatrix = current_catmatrix.catmatrix;
    all_catmatrix(1:size(current_catmatrix, 1), 1:size(current_catmatrix, 2), 1:size(current_catmatrix, 3), ss) = current_catmatrix;
    cd(prior_dir)
end

ppts_avg = nanmean(all_catmatrix, 4);
ppts_avg = nanmean(ppts_avg, 2);
ppts_avg = squeeze(ppts_avg);

%% plot region (use for V1, PEF, FEF, FFC)
ROI = [18, 198];
ROI_label1 = 'FFC_L';
ROI_label2 = 'FFC_R';

figure
subplot(2,1,1)
plot(time, ppts_avg(:, ROI(1))); 
xlim([-0.1 0.4])
title([ROI_label1 ' ROI ' num2str(ROI(1)) ' Glasser' ' : Avg across participants, across trials']);

subplot(2,1,2)
plot(time, ppts_avg(:, ROI(2))); 
xlim([-0.1 0.4])
title([ROI_label2 ' ROI ' num2str(ROI(2)) ' Glasser' ' : Avg across participants, across trials']);

%% ACC (average of 15 ROIs)
% left side
labels_L = ["L_33pr", "L_p24pr", "L_a24pr", "L_p24", "L_a24", "L_p32pr",...
          "L_a32pr", "L_d32", "L_p32", "L_s32", "L_8BM", "L_9m",...
          "L_10v", "L_10r", "L_25"];
rows_L = [];
for j = 1:length(labels_L)
    for i = 1:length(atlas.tissuelabel)
        if strcmp(atlas.tissuelabel{i}, labels_L(j))
            rows_L = [rows_L, i];
        end
    end
end

ACC_ppt_data_L = NaN(2401, length(rows_L));
for j = 1:length(rows_L)
    ACC_ppt_data_L(:,j) = ppts_avg(:,rows_L(j));
end

ACC_ppt_data_L = nanmean(ACC_ppt_data_L, 2);

% right side
labels_R = ["R_33pr", "R_p24pr", "R_a24pr", "R_p24", "R_a24", "R_p32pr",...
          "R_a32pr", "R_d32", "R_p32", "R_s32", "R_8BM", "R_9m",...
          "R_10v", "R_10r", "R_25"];
rows_R = [];
for j = 1:length(labels_R)
    for i = 1:length(atlas.tissuelabel)
        if strcmp(atlas.tissuelabel{i}, labels_R(j))
            rows_R = [rows_R, i];
        end
    end
end

ACC_ppt_data_R = NaN(2401, length(rows_R));
for j = 1:length(rows_R)
    ACC_ppt_data_R(:,j) = ppts_avg(:,rows_R(j));
end

ACC_ppt_data_R = nanmean(ACC_ppt_data_R, 2);

% plots
figure
subplot(2,1,1)
plot(time, ACC_ppt_data_L); 
xlim([-0.1 0.4])
title('ACC_L avg of 15 ROIs Glasser: Avg across participants, across trials');

subplot(2,1,2)
plot(time, ACC_ppt_data_R); 
xlim([-0.1 0.4])
title('ACC_R avg of 15 ROIs Glasser: Avg across participants, across trials');

%% PPA (average of 6 ROIs)
% left
labels_L = ["L_PHA1", "L_PHA2", "L_PHA3", "L_VMV1", "L_VMV2", "L_VMV3"];
rows_L = [];
for j = 1:length(labels_L)
    for i = 1:length(atlas.tissuelabel)
        if strcmp(atlas.tissuelabel{i}, labels_L(j))
            rows_L = [rows_L, i];
        end
    end
end

PPA_ppt_data_L = NaN(2401, length(rows_L));
for j = 1:length(rows_L)
    PPA_ppt_data_L(:,j) = ppts_avg(:,rows_L(j));
end

PPA_ppt_data_L = nanmean(PPA_ppt_data_L, 2);

% right
labels_R = ["R_PHA1", "R_PHA2", "R_PHA3", "R_VMV1", "R_VMV2", "R_VMV3"];
rows_R = [];
for j = 1:length(labels_R)
    for i = 1:length(atlas.tissuelabel)
        if strcmp(atlas.tissuelabel{i}, labels_R(j))
            rows_R = [rows_R, i];
        end
    end
end

PPA_ppt_data_R = NaN(2401, length(rows_R));
for j = 1:length(rows_R)
    PPA_ppt_data_R(:,j) = ppts_avg(:,rows_R(j));
end

PPA_ppt_data_R = nanmean(PPA_ppt_data_R, 2);

% plots
figure
subplot(2,1,1)
plot(time, PPA_ppt_data_L); 
xlim([-0.1 0.4])
title('PPA_L avg of 15 ROIs Glasser: Avg across participants, across trials');

subplot(2,1,2)
plot(time, PPA_ppt_data_R); 
xlim([-0.1 0.4])
title('PPA_R avg of 15 ROIs Glasser: Avg across participants, across trials');


%% END OF AVG ACROSS PPTS AND TRIALS 
%% fcp_3 checker: plotting relevant channels 
%  averaged across trials 

config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp4';

% check for matched MRI and MEG data
subj_match  = freeviewing_ds_pid_match(paths,step);
ssSubjPath  = @(x) paths.(subj_match.pid{x});

rangeOFsubj = 1:length(subj_match.pid);

all_ppt_data = NaN(length(subj_match.pid),1,183,2401);

for ss = rangeOFsubj % for each participant that has matched MEG/MRI data
    prior_dir = pwd();
    cd(['/home/dmatic/MEGProjects/analysis/OIRM_test_May25/analysis/' subj_match.pid{ss}])
    current_data = load('ft_meg_fullyProcessed.mat');
    current_data = current_data.data;
    ppt_data = NaN(size(current_data.trial, 2), 183, 2401);
    for i = 1:size(current_data.trial, 2)
        ppt_data(i,:,:) = current_data.trial{1}; 
    end
    ppt_data = nanmean(ppt_data, 1);
    all_ppt_data(ss,:,:,:) = ppt_data;
    time = current_data.time{1};
    cd(prior_dir)
end

all_ppt_data = nanmean(all_ppt_data, 1);
all_ppt_data = squeeze(all_ppt_data);

subplot(2,1,1)
plot(time, all_ppt_data(1, :))
xlim([-0.1 0.4])
title('V1_L (ROI 1): Plot function on average of ppt data')
subplot(2,1,2)
plot(time, all_ppt_data(181, :))
xlim([-0.1 0.4])
title('V1_R (ROI 181): Plot function on average of ppt data')
