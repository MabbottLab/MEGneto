%% Instructions
% The following script contains code used to perform statistical analaysis
% and create visual reperesnetations of results for each aim and hypothesis
% of the OIRM Faces v Scenes project

% users are required to pay attention to the paths as they must be altered
% to the meet the users' needs

%% Load config
% load the JSON file to be able to use some of the variables in analysis
% such as frequencies and times of interest for freqanalysis
config      = load_config(paths, paths.name);
config      = config.config;

%% SECTION 1: DATA NOT SEPARATED INTO CATEGORY OF VISUAL STIMULI
% this first section contains code to perform statistical analysis and
% generate figures for the first aim of the project which was to determine
% wheteher there is an increase in neural oscillatory power after
% clip/clippet changes.

%% Make some directories
% please comment out the following lines if the directories that are being
% created here already exist for you
mkdir([paths.anout_grp '/FINAL_PLOTS']);
mkdir([paths.anout_grp '/FINAL_PLOTS/groupSeparated']);
mkdir([paths.anout_grp '/FINAL_PLOTS/notgroupSeparated']);
mkdir([paths.anout_grp '/FINAL_PLOTS/groupSeparated/freqband_contribution']);
mkdir([paths.anout_grp '/FINAL_PLOTS/groupSeparated/peakActivation_times']);
mkdir([paths.anout_grp '/FINAL_PLOTS/groupSeparated/smaller_TOI_sigdiffs']);
mkdir([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output']);
mkdir([paths.anout_grp '/FINAL_PLOTS/notgroupSeparated/avgFreq_neuralResponse_perROI']);
mkdir([paths.anout_grp '/FINAL_PLOTS/notgroupSeparated/freqband_contribution']);
mkdir([paths.anout_grp '/FINAL_PLOTS/notgroupSeparated/peakActivation_times']);
mkdir([paths.anout_grp '/FINAL_PLOTS/notgroupSeparated/stats_output']);

%% Prepare input for seeing if there are neural responses to clip changes
% load the non-group separated output from fcp5
freq_structure = load([paths.anout_grp '/fcp_5_freq_structure_33msecBin.mat']);
freq_structure = freq_structure.freq_structure;

%%% configure parameters for averaging across participants in the non-visual
% content separated data. This results in a structure called freq_grandavg
% which contains a power spectrum with dims ROIs x frequency x time 
cfg = [];
cfg.keepindividual = 'no'; % no for plotting, yes for stats
cfg.toilim = 'all';
cfg.channel = 'all';
cfg.parameter = 'powspctrm';

% calculate average across particiapnts for all clips, regardless of group
% or time window of interest
freq_grandavg =    ft_freqgrandaverage(cfg, freq_structure{1},... 
                         freq_structure{2},freq_structure{3},...
                         freq_structure{4},freq_structure{5},...
                         freq_structure{6},freq_structure{7},...
                         freq_structure{8},freq_structure{9},...
                         freq_structure{10},freq_structure{11},...
                         freq_structure{12});

%%% configure parameters for averaging across particiapants in the time
% window of -1 to 0 sec (pre-stimulus). This results in a the data
% structure freq_grandavg_control which has dims ROIs x frequency x time.
% Note that the time dimensions here is smaller than in freq_grandavg as we
% are just looking for time between -1 and 0 sec.
cfg = [];
cfg.keepindividual = 'no'; % no for plotting, yes for stats
cfg.toilim = [-1 0];
cfg.channel = 'all';
cfg.parameter = 'powspctrm';

% calculate average across particiapnts for all clips, regardless of group
% but in the time window of -1 to 0 sec
freq_grandavg_control =  ft_freqgrandaverage(cfg, freq_structure{1},... 
                         freq_structure{2},freq_structure{3},...
                         freq_structure{4},freq_structure{5},...
                         freq_structure{6},freq_structure{7},...
                         freq_structure{8},freq_structure{9},...
                         freq_structure{10},freq_structure{11},...
                         freq_structure{12}); 

%%% configure parameters for averaging across particiapants in the time
% window of 0 to 1 sec (pre-stimulus). This results in a the data
% structure freq_grandavg_control which has dims ROIs x frequency x time.
% Note that the time dimensions here is smaller than in freq_grandavg as we
% are just looking for time between 0 and 1 sec.                   
cfg = [];
cfg.keepindividual = 'no'; % no for plotting, yes for stats
cfg.toilim = [0 1];
cfg.channel = 'all';
cfg.parameter = 'powspctrm';

% calculate average across particiapnts for all clips, regardless of group
% but in the time window of 0 to 1 sec
freq_grandavg_clips =    ft_freqgrandaverage(cfg, freq_structure{1},... 
                         freq_structure{2},freq_structure{3},...
                         freq_structure{4},freq_structure{5},...
                         freq_structure{6},freq_structure{7},...
                         freq_structure{8},freq_structure{9},...
                         freq_structure{10},freq_structure{11},...
                         freq_structure{12}); % average data of faces                      
                     
%% STATS ANALYSIS: for data not separated by category of visual stimuli 
% Prepare cfg to average data over time (i.e., average data in -1 to 0 sec
% and data in 0 to 1 sec over time).
cfg = [];
cfg.avgovertime = 'yes';
cfg.nanmean = 'yes';
freq_grandavg_control = ft_selectdata(cfg, freq_grandavg_control);
freq_grandavg_clips = ft_selectdata(cfg, freq_grandavg_clips);

% Something goes wrong with the dimord, fix that
freq_grandavg_control.dimord = 'subj_freq';
freq_grandavg_clips.dimord = 'subj_freq';

% Prepare cfg for stats. Indicate that we want to average across all
% frequencies and that we wish to perform a montecarlo permutation program
% followed by a one-tailed ttest. Note that the output results in one p
% value that tells us if the post-stimulus time window had increased power
% than pre-stimulus time window across all participants, frequencies, and
% ROIs
cfg = [];
cfg.avgoverfreq = 'yes';
cfg.nanmean = 'yes';
cfg.computecritval = 'yes';
cfg.computeprob = 'yes';
cfg.tail = 1;
cfg.statistic = 'depsamplesT';
cfg.correctm = 'no';
cfg.alpha = 0.05;
cfg.method = 'montecarlo';
cfg.design = [ones(1,12) ones(1,12)*2; 1:12 1:12]; 
cfg.numrandomization = 'all';
cfg.ivar = 1;
cfg.uvar = 2;
[stat] = ft_freqstatistics(cfg, freq_grandavg_clips, freq_grandavg_control); % input the post-stimulus data followed by pre-stimulus data

% Save the stats output
save([paths.anout_grp '/FINAL_PLOTS/notgroupSeparated/stats_output/avgNeuralResponse.mat'], 'stat', '-mat', '-v7.3');

%% HYPOTHESIS 1A
% plot ROI x time (avg frequency bands and avg ROIS) 

% use 'begin' and 'last' to specify which time window we are looking in (in
% this case it is -1 to 1 seconds since we care about -1 to 0 sec for
% pre-stimulus and 0 to 1 sec for post-stimulus)
begin = find(config.freqanalysis.toi == -1); % find start time index
last = find(config.freqanalysis.toi == 1.0130); % find end time index

% average data across frequencies and ROIs as that is how we performed our
% stats analysis.
twelve_ROIs_powspctrm = squeeze(nanmean(freq_grandavg.powspctrm, 2)); % average across frequencies
twelve_ROIs_powspctrm = nanmean(twelve_ROIs_powspctrm, 1); % average across ROIs
twelve_ROIs_powspctrm = twelve_ROIs_powspctrm(begin:last); % select relevant time points (-1 to 1 sec)
twelve_ROI_labels = ["V1-Left", "PEF-Left", "FEF-Left", "FFA-Left", "PPA-Left", "ACC-Left",...
          "V1-Right", "PEF-Right", "FEF-Right", "FFA-Right", "PPA-Right", "ACC-Right"]; % array of ROI labels
      
time_bin = 0.033; % indicate that data is separated by 33msec increments (this was specified in the TOI field in the JSON config)

% generate a heatmap image that will showcase power over time (-1 to 1 sec)
% and averaged across ppts, freqs, and ROIs. Then make some edits to the
% figure, add texts, change font sizes, clear it, and save it!
figure = imagesc(config.freqanalysis.toi(begin:last), min(twelve_ROIs_powspctrm):time_bin:max(twelve_ROIs_powspctrm), twelve_ROIs_powspctrm);
colormap(winter); % alter colormap 
hAxes = get(figure,'Parent');
set(hAxes, 'XTick', -1:0.2:1, 'FontSize', 20);
set(gca, 'Color', 'none');
hold on;
xline(0, 'r', 'linewidth', 5); % add vertical line at stimulus presentation
xlabel('Time (seconds)', 'FontSize', 20);
set(gca, 'YTickLabel', []);
c = colorbar;
c.Label.String = 'Power';
c.Label.FontSize = 20;
c.FontSize = 20;
title('Neural power pre-stimulus and post-stimulus', 'FontSize', 20);
saveas(gcf, [paths.anout_grp '/FINAL_PLOTS/notgroupSeparated/avgFreq_neuralResponse_33msecBin.fig']);
clf;

%% HYPOTHESIS 1B
% generate one plot per ROI - both hemispheres on one plot (code cont'd from above)
% this figure is NOT a results of stats analysis. It is simply the figure
% above, teased apart into ROIs for visualization purposes - it allows us
% to roughly see where increases in power occur for each ROI.
ROI_separated_powspctrm = squeeze(nanmean(freq_grandavg.powspctrm, 2)); % average across frequencies
time_bin = 0.033; % indicate that data is separated by 33msec increments 
twelve_ROI_labels = ["V1", "PEF", "FEF", "FFA", "PPA", "ACC",...
                    "V1", "PEF", "FEF", "FFA", "PPA", "ACC"]; % ROI labels, not hemisphere labelled as hemisphere labelling is done in the figure itself
                
% use Matlab's tiled layout funciton to plot out each ROI in desired slot
t = tiledlayout(6,2);

% general figures for each ROI and keep them on the same overall figure,
% hemisphere separated (left column for left hemisphere, right column for
% right hemisphere)
for i = 1:12
    current_data = ROI_separated_powspctrm(i,:);
    current_data = current_data(begin:last);
    
    if i == 1
        nexttile(1), imagesc(config.freqanalysis.toi(begin:last), min(current_data):time_bin:max(current_data), current_data);
        title('Left hemisphere', 'FontSize', 20);
        set(gca, 'XTickLabel', []);
    end
    if i == 2 
        nexttile(3), imagesc(config.freqanalysis.toi(begin:last), min(current_data):time_bin:max(current_data), current_data);
        set(gca, 'XTickLabel', []);
    end
    if i == 3
        nexttile(5), imagesc(config.freqanalysis.toi(begin:last), min(current_data):time_bin:max(current_data), current_data);
        set(gca, 'XTickLabel', []);
    end
    if i == 4 
        nexttile(7), imagesc(config.freqanalysis.toi(begin:last), min(current_data):time_bin:max(current_data), current_data);
        set(gca, 'XTickLabel', []);
    end
    if i == 5 
        nexttile(9), imagesc(config.freqanalysis.toi(begin:last), min(current_data):time_bin:max(current_data), current_data);
        set(gca, 'XTickLabel', []);
    end
    if i == 6 
        nexttile(11), imagesc(config.freqanalysis.toi(begin:last), min(current_data):time_bin:max(current_data), current_data);
        set(gca, 'XTick', -1:0.2:1, 'FontSize', 20);
    end
    if i == 7 
        nexttile(2), imagesc(config.freqanalysis.toi(begin:last), min(current_data):time_bin:max(current_data), current_data);
        title('Right hemisphere', 'FontSize', 20);
        set(gca, 'XTickLabel', []);
    end
    if i == 8 
        nexttile(4), imagesc(config.freqanalysis.toi(begin:last), min(current_data):time_bin:max(current_data), current_data);
        set(gca, 'XTickLabel', []);
    end
    if i == 9 
        nexttile(6), imagesc(config.freqanalysis.toi(begin:last), min(current_data):time_bin:max(current_data), current_data);
        set(gca, 'XTickLabel', []);
    end
    if i == 10 
        nexttile(8), imagesc(config.freqanalysis.toi(begin:last), min(current_data):time_bin:max(current_data), current_data);
        set(gca, 'XTickLabel', []);
    end
    if i == 11 
        nexttile(10), imagesc(config.freqanalysis.toi(begin:last), min(current_data):time_bin:max(current_data), current_data);
        set(gca, 'XTickLabel', []);
    end
    if i == 12
        % add labels and color/font edits
        nexttile(12), imagesc(config.freqanalysis.toi(begin:last), min(current_data):time_bin:max(current_data), current_data);
        set(gca, 'XTick', -1:0.2:1, 'FontSize', 20);
        title(t, 'Neural power pre-stimulus and post-stimulus', 'FontSize', 20);
        xlabel(t,'Time (seconds)', 'FontSize', 20);
        c = colorbar;
        c.Label.String = 'Power';
        c.Label.FontSize = 20;
        c.FontSize = 20;
        c.Layout.Tile = 'east';
    end
    % change the heatmap color and add a line to indicate when stimulus
    % presentation occured. 
    colormap(winter);
    hold on;
    xline(0, 'r', 'linewidth', 5); % add vertical line at stimulus presentation
    ylabel(convertStringsToChars(twelve_ROI_labels(i)), 'FontSize', 20); % place ROI label on the y-axis of each figure
    set(gca, 'YTickLabel', []); % remove y-ticks
end

% save the figure and clear uit
saveas(gcf, [paths.anout_grp '/FINAL_PLOTS/notgroupSeparated/avgFreq_neuralResponse_perROI/BothHemisphere_avgFreq_neuralResponse_33msecBin.fig']);
clf;

%% HYPOTHESIS 1C
% express power contribution from each frequency band relative to avg power
% (all clips condition) averaged across ROIs
% again, this figure is not a result of stats analysis. It is simply used
% to visually identify which frequency band contributed the most to the
% observed increases in power.

% extract average power value in the time window of interest (0 to 1s)
baseline = squeeze(nanmean(freq_grandavg_clips.powspctrm, 1)); % avg across ROIs
baseline = squeeze(nanmean(baseline, 1)); % average across freq
baseline = squeeze(nanmean(baseline)); % average across time

% get frequency data in the time window of 0-1 sec
data = squeeze(nanmean(freq_grandavg_clips.powspctrm, 1)); % avg across ROIs
data = nanmean(data, 2); % avg across time
data = (data/baseline)*100; % adjust data to reflect a percent relative to average power value

%%% First plot overall frequency (x-axis) vs power (y-axis) - across all
% ROIs
plot(config.freqanalysis.foi, data);
ax = gca;
ax.FontSize = 16;
xlabel('Frequency (Hertz)', 'FontSize', 20);
ylabel('Power increase relative to average power (%)', 'FontSize', 20);
title('Percent increase in power relative to average power for each frequency', 'FontSize', 20);
saveas(gcf, [paths.anout_grp '/FINAL_PLOTS/notgroupSeparated/freqband_contribution/percentIncrease_fromAvgVal.fig']);
clf;

%%% Now let's plot the frequency vs power separated by each ROI to identify
% if there is any difference in freq band contribution between ROIs
data = squeeze(nanmean(freq_grandavg_clips.powspctrm, 3)); % avg across time
twelve_ROIs_labels = ["V1-Left", "PEF-Left", "FEF-Left", "FFA-Left", "PPA-Left", "ACC-Left",...
                      "V1-Right", "PEF-Right", "FEF-Right", "FFA-Right", "PPA-Right", "ACC-Right"];

% use MATLAB's tiled layout function to facilitate palcement of plots 
figure;
t = tiledlayout(6,2);

for i = 1:12
    % extract the data you need
    data_subset = data(i, :);
    data_subset = (data_subset/baseline)*100; % specify % increase
    
    % plot tiles
    nexttile, plot(config.freqanalysis.foi, data_subset);
    title(twelve_ROIs_labels(i), 'FontSize', 20);
    ax = gca;
    ax.FontSize = 12;
    
    if i == 12
        title(t, 'Percent increase in power relative to average power for each frequency', 'FontSize', 20);
        ylabel(t,'Power increase relative to average power (%)', 'FontSize', 20);
        xlabel(t,'Frequency (Hertz)', 'FontSize', 20);
    end
    
end

% save and clear figure
saveas(gcf, [paths.anout_grp '/FINAL_PLOTS/notgroupSeparated/freqband_contribution/ROI_separated_percentIncrease_fromAvgVal.fig']);
clf;

%% SECTION 2: DATA SEPARATED INTO CATEGORY OF VISUAL STIMULI (FACES V OTHER)
% this second section contains code to perform statistical analysis and
% generate figures for the second aim of the project which was to determine
% whether visual content containing human faces had an increase in power
% across all frequency bands and ROIs (save the PPA) in a given time
% window.

% The initial hypothesis was that there was an increase in power across all
% freq band and ROIs except the PPA in the averaged time window of 0-1
% seconds. Then, exploratory analysis was done to investigate this for
% smaller time bins (increments of 50msec).

%% Prepare input for visual stimuli-separated analysis
% faces data
freq_structure_faces = load([paths.anout_grp '/fcp_5_freq_structure_faces.mat']);
% other visual content data
freq_structure_other = load([paths.anout_grp '/fcp_5_freq_structure_other.mat']);
freq_structure_faces = freq_structure_faces.freq_structure;
freq_structure_other = freq_structure_other.freq_structure;

% configure some settings for selecting data
cfg = [];
cfg.keepindividual = 'yes'; % CHANGE THIS FIELD DEPENDING ON IF DOING STATS ('yes') OR IF DOING PLOTS ('no')
cfg.foilim = 'all';
cfg.toilim = 'all';
cfg.channel = 'all';
cfg.parameter = 'powspctrm';

% average across all ppts for faces data
freq_grandavg_faces =    ft_freqgrandaverage(cfg, freq_structure_faces{1},... 
                         freq_structure_faces{2},...
                         freq_structure_faces{3}, freq_structure_faces{4},...
                         freq_structure_faces{5}, freq_structure_faces{6},...
                         freq_structure_faces{7}, freq_structure_faces{8},...
                         freq_structure_faces{9}, freq_structure_faces{10},...
                         freq_structure_faces{11}, freq_structure_faces{12});

% average across all ppts for other visual content
freq_grandavg_other =    ft_freqgrandaverage(cfg, freq_structure_other{1},...
                         freq_structure_other{2},...
                         freq_structure_other{3}, freq_structure_other{4},...
                         freq_structure_other{5}, freq_structure_other{6},...
                         freq_structure_other{7}, freq_structure_other{8},...
                         freq_structure_other{9}, freq_structure_other{10},...
                         freq_structure_other{11}, freq_structure_other{12});                    
                                         
%% STATS ANALYSIS: for data separated into category of visual stimuli
%% Stats for data averaged across 0 to 1 seconds (initial hypothesis)
% create labels for each frequency band of interest
freq_labels = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"];
freq_ranges = {[0 4], [4 8], [8 12], [13 30], [30 60], [60 100]};
all_freq_ranges = cell(2,6);

% perform statistics for each frequency band and each ROI, averaged across 
% 0-1 second and averaged across ppts. This analysis reveals where faces 
% had increased power compared to other visual content

% Did NOT correct for multiple comparisons since a priori hypotheses were
% established for each frequency band and ROI
for i = 1:length(freq_labels)
    cfg = [];
    cfg.frequency = freq_ranges{i};
    cfg.avgoverfreq = 'yes';
    cfg.nanmean = 'yes';
    cfg.latency = [0 1];
    cfg.avgovertime = 'yes';
    cfg.computecritval = 'yes';
    cfg.computeprob = 'yes';
    cfg.tail = 1;
    cfg.statistic = 'depsamplesT';
    cfg.correctm = 'no';
    cfg.alpha = 0.05;
    cfg.method = 'montecarlo';
    cfg.design = [ones(1,12) ones(1,12)*2; 1:12 1:12]; 
    cfg.numrandomization = 'all';
    cfg.ivar = 1;
    cfg.uvar =2;
    [stat] = ft_freqstatistics(cfg, freq_grandavg_faces, freq_grandavg_other);
    
    % discount initially NaN values from significance vectors
    stat.prob(find(isnan(stat.stat))) = NaN;
    stat.prob(find(isnan(stat.stat))) = NaN;
    stat.mask(find(isnan(stat.stat))) = 0;
    stat.mask(find(isnan(stat.stat))) = 0;
    
    % store output
    all_freq_ranges{1, i} = freq_labels(i);
    all_freq_ranges{2, i} = stat;
end

% add in significance vector (p<0.05) to mark which freq bands and ROIs
% have significant differences
for i = 1:length(all_freq_ranges)
    ROI_splitter = cell(12,2);
    ROI_splitter(:,1) = all_freq_ranges{2,i}.label;
    
    for j = 1:size(all_freq_ranges{2,i}.prob, 1)
        significant_indices = find(all_freq_ranges{2,i}.prob(j,:,:) <= 0.05);
        if ~isempty(significant_indices)
            ROI_splitter(j, 2) = {significant_indices};
        end 
    end
   
    all_freq_ranges{3,i} = ROI_splitter;
    
    if isempty(all_freq_ranges{3,i})
        all_freq_ranges{3,i} = NaN;
    end
end

% save output
save([paths.anout_grp '/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_0_to_1.mat'], 'all_freq_ranges', '-mat', '-v7.3');

%% STATS ANALYSIS: for data separated into category of visual stimuli
%% Stats analysis for exploratory analysis - time bins of 50msec
% create labels for each frequency band of interest
freq_labels = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"];
freq_ranges = {[0 4], [4 8], [8 12], [13 30], [30 60], [60 100]};
all_freq_ranges = cell(2,6);

% perform statistics for each frequency band and ROI, averaged across a 
% user-specified 50msec time bin (see line of code with comment 
% 'CHANGE THIS FIELD'), averaged across ppts.This analysis reveals where 
% faces had significantly increased power compared to other visual content.

% Did NOT correct for multiple comparisons since this is explatory.

for i = 1:length(freq_labels)
    cfg = [];
    cfg.frequency = freq_ranges{i};
    cfg.avgoverfreq = 'yes';
    cfg.nanmean = 'yes';
    cfg.latency = [0.9 0.95]; % CHANGE THIS FIELD to the 50msec timewindow of interest!!! 
    cfg.avgovertime = 'yes';
    cfg.computecritval = 'yes';
    cfg.computeprob = 'yes';
    cfg.tail = 1;
    cfg.statistic = 'depsamplesT';
    cfg.correctm = 'no';
    cfg.alpha = 0.05;
    cfg.method = 'montecarlo';
    cfg.design = [ones(1,12) ones(1,12)*2; 1:12 1:12]; 
    cfg.numrandomization = 'all';
    cfg.ivar = 1;
    cfg.uvar =2;
    [stat] = ft_freqstatistics(cfg, freq_grandavg_faces, freq_grandavg_other);
   
    % discount initially NaN values from significance vectors
    stat.prob(find(isnan(stat.stat))) = NaN;
    stat.prob(find(isnan(stat.stat))) = NaN;
    stat.mask(find(isnan(stat.stat))) = 0;
    stat.mask(find(isnan(stat.stat))) = 0;
    
    % store output
    all_freq_ranges{1, i} = freq_labels(i);
    all_freq_ranges{2, i} = stat;
end

% add in significance vector
for i = 1:length(all_freq_ranges)
    ROI_splitter = cell(12,2);
    ROI_splitter(:,1) = all_freq_ranges{2,i}.label;
    
    for j = 1:size(all_freq_ranges{2,i}.prob, 1)
        significant_indices = find(all_freq_ranges{2,i}.prob(j,:,:) <= 0.05);
        if ~isempty(significant_indices)
            ROI_splitter(j, 2) = {significant_indices};
        end 
    end
   
    all_freq_ranges{3,i} = ROI_splitter;
    
    if isempty(all_freq_ranges{3,i})
        all_freq_ranges{3,i} = NaN;
    end
end

% save output
save([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_900_to_950msec.mat'], 'all_freq_ranges', '-mat', '-v7.3'); % CHANGE THIS FIELD to the output name that matches your timewindow of interest!!! 
%% Load stats output from the 50msec time bin increments 
% Load all outputs from the stats analysis above. The commented ones are
% there as an example of the time bins that were previously done in the
% OIRm analysis. You may wish to use different ones.

all_freq_ranges_Montecarlo_onetail_0to1 = load([paths.anout_grp '/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_0_to_1.mat']);
all_freq_ranges_Montecarlo_onetail_0to50msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_0_to_50msec.mat']);
all_freq_ranges_Montecarlo_onetail_50to100msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_50_to_100msec.mat']);
all_freq_ranges_Montecarlo_onetail_100to150msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_100_to_150msec.mat']);
all_freq_ranges_Montecarlo_onetail_150to200msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_150_to_200msec.mat']);
all_freq_ranges_Montecarlo_onetail_200to250msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_200_to_250msec.mat']);
all_freq_ranges_Montecarlo_onetail_250to300msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_250_to_300msec.mat']);
all_freq_ranges_Montecarlo_onetail_300to350msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_300_to_350msec.mat']);
all_freq_ranges_Montecarlo_onetail_350to400msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_350_to_400msec.mat']);
all_freq_ranges_Montecarlo_onetail_400to450msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_400_to_450msec.mat']);
all_freq_ranges_Montecarlo_onetail_450to500msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_450_to_500msec.mat']);
all_freq_ranges_Montecarlo_onetail_500to550msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_500_to_550msec.mat']);
all_freq_ranges_Montecarlo_onetail_550to600msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_550_to_600msec.mat']);
all_freq_ranges_Montecarlo_onetail_600to650msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_600_to_650msec.mat']);
all_freq_ranges_Montecarlo_onetail_650to700msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_650_to_700msec.mat']);
all_freq_ranges_Montecarlo_onetail_700to750msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_700_to_750msec.mat']);
all_freq_ranges_Montecarlo_onetail_750to800msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_750_to_800msec.mat']);
all_freq_ranges_Montecarlo_onetail_800to850msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_800_to_850msec.mat']);
all_freq_ranges_Montecarlo_onetail_850to900msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_850_to_900msec.mat']);
all_freq_ranges_Montecarlo_onetail_900to950msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_900_to_950msec.mat']);
all_freq_ranges_Montecarlo_onetail_950to1000msec = load([paths.anout_grp '/FINAL_PLOTS/groupSeparated/stats_output/fcp_5_5_MonteCarlo_all_freq_ranges_onetail_uncorrected_avg_over_950_to_1000msec.mat']);

%% HYPOTHESIS 2: Faces v Other differences in power
% This code is used to generate a plot that show cases where there is a
% significant difference in power in a given ROI and freq band for faces
% compared to other visual content. The user must changed the time window
% and the structure being analysed (of the ones loaded in the prior
% section). To identify lines of code that must be changed by the user,
% look for comments that say 'CHAGNE THIS FIELD'. Also, note that the
% smaller time bin (50msec) separated figure was generated in an excel
% sheet that is located on the Mabbott Lab Sharepoint.

% specify labels for frequencies and ROIs
freq_labels = ["delta", "theta", "alpha", "beta", "low-gamma", "high-gamma"];
labels = ["V1-Left", "PEF-Left", "FEF-Left", "FFA-Left", "PPA-Left", "ACC-Left",...
          "V1-Right", "PEF-Right", "FEF-Right", "FFA-Right", "PPA-Right", "ACC-Right"];

% specify which stat results are being used      
all_freq_ranges = all_freq_ranges_Montecarlo_onetail_0to1; % CHANGE THIS FIELD depending on time window of interest

% specify time window of interest
time_window = "0-to-1sec"; % CHANGE THIS FIELD depending on time window of interest

% create empty structures to hold output data
all_data = zeros(6,12);
all_mask = zeros(6,12);

for i = 1:length(all_freq_ranges.all_freq_ranges(3,:))
    mask = squeeze(all_freq_ranges.all_freq_ranges{2,i}.mask);
    data = squeeze(all_freq_ranges.all_freq_ranges{2,i}.prob);
    all_data(i,:) = data.';
    all_mask(i,:) = mask;
end

% generate figure and mark areas where significant differences occured with
% a red star
figure = imagesc(1:12, 1:6, all_data);
colormap(flipud(parula));
hold on;
for i = 1:height(all_mask)
    for j = 1:length(all_mask)
        if all_mask(i,j) == 1
            plot(j,i,'*r');
            hold on;
        end
    end
end

% add labels and change font
hAxes = get(figure,'Parent');
set(hAxes, 'YTick', 1:6, 'YTickLabel', freq_labels, 'FontSize', 12)
set(hAxes, 'XTick', 1:12, 'XTickLabel', labels, 'FontSize', 12)
title(["Difference in power for faces vs. other visual stimulus " time_window], 'FontSize', 12);
c = colorbar;
c.Label.String = 'p-value';
c.Label.FontSize = 12;

% save and clear figure
saveas(gcf, [paths.anout_grp '/FINAL_PLOTS/groupSeparated/smaller_TOI_sigdiffs/' convertStringsToChars(time_window) '.fig'], 'fig');
clf;

%% UNUSED CODE
% This code was not used but remains here in the case that it is useful in
% the future. The code can be used to identify the highest power value in a
% given region of interest across frequency bands. 

%% Peak activation times for data not separated into category of visual stimuli
% twelve_ROIs_powspctrm = squeeze(nanmean(freq_grandavg.powspctrm, 2)); % average across frequencies
% % select times!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% begin = find(config.freqanalysis.toi == 0); % find start time index
% last = find(config.freqanalysis.toi == 1); % find end time index
% 
% all_max = array2table(zeros(12,2), 'VariableNames', {'X', 'Y'});
% all_max.Properties.RowNames = twelve_ROI_labels;
% for k = 1:12 % number of ROIs
%     % extract peak coordinates
%     [x_all, y_all] = max(twelve_ROIs_powspctrm(k, begin:last)); % this data is already 12 x 81 (ROI x time)
%     
%     % store in table
%     time = config.freqanalysis.toi(begin:last);
%     all_max(k,1) = {time(y_all)};
%     all_max(k,2) = {x_all};
% end
% % sort rows
% all_max = sortrows(all_max,1);
% 
% save(paths.anout_grp '/FINAL_PLOTS/notgroupSeparated/peakActivation_times/peaktimes.mat'], 'all_max', '-mat', '-v7.3');
% 
% % create plot of peak activation times
% order = [];
% for i = 1:12
%     index = find(all_max.Properties.RowNames(i) == twelve_ROI_labels);
%     order = [order index];
% end
% 
% for i = 1:12
%     subplot(12,1,i);
%     % plot waveform and mark the peak position
%     figure = plot(config.freqanalysis.toi(begin:last), twelve_ROIs_powspctrm(order(i),begin:last));
%     hold on;
%     plot(all_max.X(i), all_max.Y(i), '*r');
%     hold on;
%     text(all_max.X(i) + 0.03, all_max.Y(i) - 0.003, [num2str(all_max.X(i)) ' sec'], 'FontSize', 12);
%     hold on;
%     
%     % set and adjust y-label
%     ylabel(all_max.Properties.RowNames{i});
%     set(gca, 'YTickLabel', []);
%     hYLabel = get(gca,'YLabel');
%     set(hYLabel,'rotation',0, 'VerticalAlignment', 'middle', 'FontSize', 12)
%     hYLabel.Position(1) = config.freqanalysis.toi(begin) - 0.03;
%     hXLabel = get(gca, 'XLabel');
%     set(hXLabel, 'FontSize', 12);
%     
%     % set x-label
%     if i == 12
%         xlabel('Time (seconds)', 'FontSize', 12);
%     end
%     
%     % supress x-ticks for all subplots but the last one
%     if i ~= 12
%         set(gca, 'XTickLabel', [])
%     end
%     
%     % set title at the top subplot
%     if i == 1
%         title('Peak activation times of ROIs across all clips', 'FontSize', 12);
%     end
% end
% 
% saveas(gcf, [paths.anout_grp '/FINAL_PLOTS/notgroupSeparated/peaktimes.fig']);
% clf;

%% Peak activation times for data separated into category of visual stimuli in each ROI (avg across frequency bands)
% freq_labels = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"];
% freq_ranges = {[3 4], [4 8], [8 12], [13 30], [30 60], [60 100]};
% twelve_ROIs_labels = ["V1-Left", "PEF-Left", "FEF-Left", "FFA-Left", "PPA-Left", "ACC-Left",...
%                       "V1-Right", "PEF-Right", "FEF-Right", "FFA-Right", "PPA-Right", "ACC-Right"];
% 
% % create cell to store data
% facesdata_powspctrm = {};
% otherdata_powspctrm = {};
% 
% begin = find(config.freqanalysis.toi == 0); % find start time index
% last = find(config.freqanalysis.toi == 1); % find end time index
% 
% for i = 1:length(freq_ranges) % number of freq bands
%     faces_max = array2table(zeros(12,2), 'VariableNames', {'X', 'Y'});
%     other_max = array2table(zeros(12,2), 'VariableNames', {'X', 'Y'});
%     faces_max.Properties.RowNames = twelve_ROIs_labels;
%     other_max.Properties.RowNames = twelve_ROIs_labels;
%     
%     cfg = [];
%     cfg.avgoverfreq = 'yes';
%     cfg.nanmean = 'yes';
%     data1 = ft_selectdata(cfg, freq_grandavg_faces);
%     data2 = ft_selectdata(cfg, freq_grandavg_other);
%     data1_powspctrm = squeeze(data1.powspctrm);
%     data2_powspctrm = squeeze(data2.powspctrm);
%     
%     for k = 1:12 % number of ROIs
%         % extract peak coordinates
%         [x_face, y_face] = max(data1_powspctrm(k, begin:last));
%         [x_other, y_other] = max(data2_powspctrm(k, begin:last));
%         
%         % store in table
%         time = config.freqanalysis.toi(begin:last);
%         faces_max(k,1) = {time(y_face)};
%         faces_max(k,2) = {x_face};
%         other_max(k,1) = {time(y_other)};
%         other_max(k,2) = {x_other};
%     end
%     % sort rows
%     faces_max = sortrows(faces_max,1);
%     other_max = sortrows(other_max,1);
% end
% 
% save([paths.anout '/group/FINAL_PLOTS/groupSeparated/peakActivation_times/twelveROIs_faces_peaktime.mat'], 'faces_max', '-mat', '-v7.3');
% save([paths.anout '/group/FINAL_PLOTS/groupSeparated/peakActivation_times/twelveROIs_other_peaktime.mat'], 'other_max', '-mat', '-v7.3');
% 
% % avg data over frequency
% faces = squeeze(nanmean(freq_grandavg_faces.powspctrm, 2));
% other = squeeze(nanmean(freq_grandavg_other.powspctrm, 2));
% 
% % create plot for FACES
% faces_order = [];
% for i = 1:12
%     index = find(faces_max.Properties.RowNames(i) == twelve_ROIs_labels);
%     faces_order = [faces_order index];
% end
% 
% for i = 1:12
%     subplot(12,1,i);
%     % plot waveform and mark the peak position
%     plot(time, faces(faces_order(i),begin:last));
%     hold on;
%     plot(faces_max.X(i), faces_max.Y(i), '*r');
%     hold on;
%     text(faces_max.X(i) + 0.03, faces_max.Y(i) - 0.01, [num2str(faces_max.X(i)) ' sec'], 'FontSize', 12);
%     hold on;
%     
%     % set and adjust y-label
%     ylabel(faces_max.Properties.RowNames{i});
%     set(gca, 'YTickLabel', []);
%     hYLabel = get(gca,'YLabel');
%     set(hYLabel,'rotation',0, 'VerticalAlignment', 'middle', 'FontSize', 12)
%     hYLabel.Position(1) = config.freqanalysis.toi(begin)-0.02;
%     hXLabel = get(gca, 'XLabel');
%     set(hXLabel, 'FontSize', 12);
%     
%     % supress x-ticks for all subplots but the last one
%     if i ~= 12
%         set(gca, 'XTickLabel', [])
%     end
%     
%     % set title at the top subplot
%     if i == 1
%         title('Peak activation times of ROIs for clips containing faces', 'FontSize', 12);
%     end
%     
%     if i == 12
%         xlabel('Time (seconds)', 'FontSize', 12);
%     end
% end
% saveas(gcf, [paths.anout_grp '/FINAL_PLOTS/groupSeparated/peakActivation_times/faces_peakTimes.fig']);
% clf;
% 
% % create plot for OTHER
% other_order = [];
% for i = 1:12
%     index = find(other_max.Properties.RowNames(i) == twelve_ROIs_labels);
%     other_order = [other_order index];
% end
% 
% for i = 1:12
%     subplot(12,1,i);
%     % plot waveform and mark the peak position
%     plot(time, other(other_order(i), begin:last));
%     hold on;
%     plot(other_max.X(i), other_max.Y(i), '*r');
%     hold on;
%     text(other_max.X(i) + 0.03, other_max.Y(i) - 0.01, [num2str(other_max.X(i)) ' sec'], 'FontSize', 12);
%     hold on;
%     
%     % set and adjust y-label
%     ylabel(other_max.Properties.RowNames{i});
%     set(gca, 'YTickLabel', []);
%     hYLabel = get(gca,'YLabel');
%     set(hYLabel,'rotation',0, 'VerticalAlignment', 'middle', 'FontSize', 12)
%     hYLabel.Position(1) = config.freqanalysis.toi(begin)-0.02;
%     hXLabel = get(gca, 'XLabel');
%     set(hXLabel, 'FontSize', 12);
%     
%     % supress x-ticks for all subplots but the last one
%     if i ~= 12
%         set(gca, 'XTickLabel', [])
%     end
%     
%     % set title at the top subplot
%     if i == 1
%         title('Peak activation times of ROIs for clips containing no faces', 'FontSize', 12);
%     end
%     
%     % x-axis label
%     if i == 12
%         xlabel('Time (seconds)', 'FontSize', 12);
%     end
% end
% saveas(gcf, [paths.anout_grp '/FINAL_PLOTS/groupSeparated/peakActivation_times/other_peakTimes.fig']);
% clf;



