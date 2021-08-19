%% SCRAP CODE
% This script contains several scrap bits of code that may come in handy
% and are left here for future convenience. This code was not used for
% analysis.

%% Setting up paths
addpath('/mnt/hpc-megneto/fieldtrip')
ft_defaults
addpath(genpath('/mnt/hpc-megneto/MEGneto'))
dunja_MEG_code = '/home/dmatic/MEGcode';
addpath(dunja_MEG_code)

%% Loading data and using ft_timelockgrandaverage
%% Run 1
config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp3';

% check for matched MRI and MEG data
subj_match  = freeviewing_ds_pid_match(paths,step);
ssSubjPath  = @(x) paths.(subj_match.pid{x});

rangeOFsubj = 1:length(subj_match.pid);

%all_ppt_data = NaN(length(subj_match.pid),1,183,2401);

tl_output = cell(1,height(subj_match));

for ss = rangeOFsubj % for each participant that has matched MEG/MRI data
    prior_dir = pwd();
    cd(['/home/dmatic/MEGProjects/analysis/OIRMProcessingRun1_Jun16/analysis/' subj_match.pid{ss}])
    
    current_data = load('ft_meg_fullyProcessed.mat');
    current_data = current_data.data;
    
    cfg = [];
    cfg.channel = 'MEG';
    tl = ft_timelockanalysis(cfg, current_data);
    tl_output{ss} = tl;
    
    cd(prior_dir)
end

%% Run 2
config2      = load_config(paths, paths.name);
config2      = config2.config;
step        = 'fcp3';

% check for matched MRI and MEG data
subj_match2  = freeviewing_ds_pid_match(paths,step);
ssSubjPath2  = @(x) paths.(subj_match2.pid{x});
subj_match2.ds(~cellfun(@isempty,regexp(subj_match2.ds,'Free.ds','match')))=[];
% subj_match2(1,:) = [];
% subj_match2(2,:) = [];
% subj_match2(3,:) = [];
% subj_match2(4,:) = [];
% subj_match2(5,:) = [];
% subj_match2(6,:) = [];
% subj_match2(7,:) = [];
% subj_match2(8,:) = [];
% subj_match2(11,:) = [];

rangeOFsubj2 = 1:length(subj_match2.pid);

%all_ppt_data = NaN(length(subj_match.pid),1,183,2401);

tl_output2 = cell(1,height(subj_match2));

for ss = rangeOFsubj2 % for each participant that has matched MEG/MRI data
    prior_dir = pwd();
    cd(['/home/dmatic/MEGProjects/analysis/OIRM_test_May25/analysis/' subj_match2.pid{ss}])
    
    current_data2 = load('ft_meg_fullyProcessed.mat');
    current_data2 = current_data2.data;
    
    cfg = [];
    cfg.channel = 'MEG';
    tl2 = ft_timelockanalysis(cfg, current_data2);
    tl_output2{ss} = tl2;
    
    cd(prior_dir)
end
%% Grand average and plots
% take average across many participants
cfg =[];
cfg.channel = 'MEG';
tl_grandaverage = ft_timelockgrandaverage(cfg, tl_output{1},tl_output{2}, tl_output{3}, tl_output{4}, tl_output{5},...
                                          tl_output{6}, tl_output{7}, tl_output{8}, tl_output{9}, tl_output{10},...
                                          tl_output2{1}, tl_output2{2}, tl_output2{3}, tl_output2{4}, tl_output2{5},...
                                          tl_output2{6}, tl_output2{7}, tl_output2{8}, tl_output2{9}, tl_output2{10},...
                                          tl_output2{11});

% sensor space paper 
cfg = [];
cfg.xlim = [-0.1, 0.4];
cfg.channel = {'MLP11', 'MLP12', 'MLP21', 'MLP22', 'MLP31', 'MLP32',...
    'MLO11', 'MLO12', 'MLO21', 'MLO22', 'MLO31', 'MLO32', 'MRP11', 'MRP12',...
    'MRP21', 'MRP22', 'MRP31', 'MRP32', 'MRO11', 'MRO12', 'MRO21', 'MRO22',...
    'MRO31', 'MRO32', 'MZP01', 'MZP02', 'MZO01', 'MZO02'};
ft_singleplotER(cfg, tl_grandaverage)
title('parietooccipital channels: ft_singleplotER on tl_grandaverage');

                                      
% fusiform gyrus is occipital and temporal lobe
% right
cfg = [];
cfg.xlim = [-0.1, 0.4];
cfg.channel = {'MRO*', 'MRT1*', 'MRT2*', 'MRT3*', 'MRT42', 'MRT43', 'MRT44'};
ft_singleplotER(cfg, tl_grandaverage)
title('RIGHT occipito-temporal channels (fusiform gyrus?): ft_singleplotER on tl_grandaverage');
% left
cfg = [];
cfg.xlim = [-0.1, 0.4];
cfg.channel = {'MLO*', 'MLT*'};
ft_singleplotER(cfg, tl_grandaverage)
title('LEFT occipito-temporal channels (fusiform gyrus?): ft_singleplotER on tl_grandaverage');
% all
cfg = [];
cfg.xlim = [-0.1, 0.4];
cfg.channel = {'MRO*', 'MRT1*', 'MRT2*', 'MRT3*', 'MRT42', 'MRT43', 'MRT44', 'MLO*', 'MLT*'};
ft_singleplotER(cfg, tl_grandaverage)
title('ALL occipito-temporal channels (fusiform gyrus?): ft_singleplotER on tl_grandaverage');
                                      
% look at select channels in single plots
% Parietal
cfg = [];
cfg.xlim = [-0.1, 0.4];
cfg.channel = {'MLP*', 'MZP*', 'MRP*'};
ft_singleplotER(cfg, tl_grandaverage)
title('Parietal channels: ft_singleplotER on tl_grandaverage');

% Temporal
cfg = [];
cfg.xlim = [-0.1, 0.4];
cfg.channel = {'MLT*', 'MRT1*', 'MRT2*', 'MRT3*', 'MRT42', 'MRT43', 'MRT44'};
ft_singleplotER(cfg, tl_grandaverage)
title('Temporal channels: ft_singleplotER on tl_grandaverage');

% Frontal
cfg = [];
cfg.xlim = [-0.1, 0.4];
cfg.channel = {'MLF*', 'MZF*', 'MRF*'};
ft_singleplotER(cfg, tl_grandaverage)
title('Frontal channels: ft_singleplotER on tl_grandaverage');

% Occipital
cfg = [];
cfg.xlim = [-0.1, 0.4];
cfg.channel = {'MLO*', 'MZO*', 'MRO*'};
ft_singleplotER(cfg, tl_grandaverage)
title('Occipital channels: ft_singleplotER on tl_grandaverage');

% Central
cfg = [];
cfg.xlim = [-0.1, 0.4];
cfg.channel = {'MLC*', 'MZC*', 'MRC*'};
ft_singleplotER(cfg, tl_grandaverage)
title('Central channels: ft_singleplotER on tl_grandaverage');

% plot all channels and topography 
cfg = [];
cfg.showlabels = 'yes';
cfg.showoutline = 'yes';
cfg.xlim = [-0.1, 0.4];
cfg.layout = 'CTF151_helmet.mat';
ft_multiplotER(cfg, tl_grandaverage);
% ft_topoplotER(cfg, tl_grandaverage);

% plot individual ppt activity too
figure; 
for isub = 1:length(subj_match.pid)
    prior_dir = pwd();
    cd(['/home/dmatic/MEGProjects/analysis/OIRM_test_May25/analysis/' subj_match.pid{isub}])
    current_data = load('ft_meg_fullyProcessed.mat');
    current_data = current_data.data; 
    cd(prior_dir);
    
    cfg = [];
    cfg.xlim = [-0.1, 0.4];
    cfg.channel = {'MLO*', 'MZO*', 'MRO*'};
    
    subplot(3,4,1); ft_singleplotER(cfg, current_data);
end