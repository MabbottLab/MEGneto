function fcp_4_beamforming(paths)

% FCP_4_BEAMFORMING carries out beamforming and source projection on
% cleaned data. 
% 
% NOTES:
%   - Ensure that subj_fcp4.csv is populated with the subject IDs of
%   participants you want to include after checking over initial results. 
%   - Caution: unit conversion is sensitive in this step b/w cm and mm. 
%   - Assumes that anatomical MRI is already aligned with MEG data. If not
%   aligned, use ft_volumerealign function. 
%
% INPUTS:
%   paths               =   struct defining paths to data, participant
%                           folders, analysis folders, config files, etc. 
%
% OUTPUTS:
%   .mat
%       - catmatrix:    num_samples x num_trials x num_nodes
%       - srate:        sampling rate
%       - coords:       coordinates of source points w/in
%
% See also: 
%
% Last updated by: Julie Tseng, 2020-01-08
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.

%% SET UP LOGGING FILE

right_now = clock;
log_filename = [paths.conf_dir '/log_' sprintf('%d%d%d', right_now(1:3))];
diary(log_filename)

fprintf('\n\n%d:%d:%02.f       Now running **%s**.\n', ...
    right_now(4:6), mfilename)

%% SETUP

% load config JSON with analysis parameters
config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp4';

% check for matched MRI and MEG data
subj_match  = ds_pid_match(paths,step);
ssSubjPath  = @(x) paths.(subj_match.pid{x});
if isempty(subj_match) % if there are no full sets of data
    error('No participants selected')
end

% save the final participants list to the subj_match_fcp1 CSV
write_match_if_not_empty(paths,step);

% check for multiple *.ds folders for participants
if length(unique(subj_match.pid)) ~= length(subj_match.pid)
    error('More than one ds per participant!')
end
% save the final participants list to the subj_match_fcp1 CSV
write_match_if_not_empty(paths,step);

% check for multiple *.ds folders for participants
if length(unique(subj_match.pid)) ~= length(subj_match.pid)
    error('More than one ds per participant!')
end

% initialize output files
    % images
%         fcp4_output.fig_headsurf  = 'headsurf.png';

%% PARTICIPANT MODELS

rangeOFsubj = 1:length(subj_match.ds);

for ss = rangeOFsubj
    %%% FOR EACH PARTICIPANT --------------------------------------------------
    fprintf('\nWorking on subject %s! \n', subj_match.pid{ss});
    
    %%% Load Headmodel and surface 8k -------------------------------
    load(fullfile(ssSubjPath(ss),sprintf('/%s_sourcemodel_4k.mat',subj_match.pid{ss})),'-mat', 'sourcemodel_4k');
    load(fullfile(ssSubjPath(ss),sprintf('/%s_hdm.mat',subj_match.pid{ss})), '-mat','headmodel');
 
     %%% LOAD MEG DATA ---------------------------------------------------------
    load([ssSubjPath(ss) '/ft_meg_fullyProcessed.mat'],'-mat','data');

    %% VISUALIZATION: Headmodel and Surfaces --------------------------------
%     close all;    
%     figure;
%         % make the headmodel surface transparent
%         ft_plot_headmodel(headmodel, 'edgecolor', 'none'); alpha 0.5
%         ft_plot_mesh(ft_convert_units(sourcemodel_8k, 'cm'),'vertexcolor',sourcemodel_8k.sulc);
%         ft_plot_sens(data.grad);
%         view([0 -90 0])
%                 
%     saveas(gcf,fullfile(ssSubjPath(ss),'/surf_head.png'))
%     savefig(fullfile(ssSubjPath(ss),'/surf_head.fig'))
    
    sourcemodel_4k = ft_convert_units(sourcemodel_4k, 'cm');
    %%% Compute the leadfield ----------------------------------------------------
    cfg             = [];
    cfg.grid        = sourcemodel_4k;
    cfg.headmodel   = headmodel;
    cfg.channel     = {'MEG'};
    lf              = ft_prepare_leadfield(cfg, data);
    
    meg = ft_channelselection({'MEG'},data);
%%
    %%% Spectral analysis ---------------------------------------------------------
    % Compute sensor level Fourier spectra, to be used for cross-spectral density computation
    cfg             = [];
    cfg.method      = 'mtmfft';%config.beamforming.freqDomain.method; %'mtmfft'; % analysis an entire spectrum for entire data length (mulitaper freq transformation)
    cfg.output      = 'fourier'; %config.beamforming.freqDomain.output; %'fourier'; % return complex Fourier-spcetra ; powandcsd return the power and cross-spectra
%     cfg.taper       = 'dpss'; %config.beamforming.freqDomain.taper; %'dpss'; %default 
    cfg.channel     = meg;
    cfg.tapsmofrq   = 9; %config.beamforming.freqDomain.tapsmofrq; %1;
    cfg.foi         = 39; %config.beamforming.freqDomain.foi; %10 alpha ; 5 theta 
    cfg.keeptrials  = 'yes'; %config.beamforming.options.keeptrials; %'yes'
    %wpli only
    % cfg.t_ftimwin = ones(length(cfg.foi),1).*0.5; %length of time window in seconds = 0.5s
    % cfg.toi     = 0:0.05:1.5; %the times of which the analysis windows should be centered (in seconds) from 0 to 6.99s in steps of 0.05   
    freq           = ft_freqanalysis(cfg,data);



    %%% Forward Model ---------------------------------------------------------
    % Use a cortical sheeet based source model, in which individual dipole locations 
    % are constrated to cortical sheet. 

    % Source reconstructrion    
    cfg                 = [];
    cfg.frequency       = freq.freq;
    cfg.method           = 'pcc'; %config.beamforming.method; % eloreta or pcc (coh)
    cfg.grid            = lf;
    %cfg.sourcemodel      =lf; 
    cfg.headmodel       = headmodel;
    cfg.keeptrials      = 'yes';
    cfg.pcc.lambda      = '10%';
    cfg.pcc.projectnoise = 'yes';
    cfg.pcc.fixedori    ='yes';  
    source = ft_sourceanalysis(cfg, freq);


    % save results
    save(fullfile(ssSubjPath(ss),sprintf('/%s_meg_sourcemodel_4k_%s',subj_match.pid{ss},'gamma1')), 'source');
    save(fullfile(ssSubjPath(ss),sprintf('/%s_meg_freq_%s',subj_match.pid{ss},'gamma1')), 'freq');
%     save([ssSubjPath(ss) '/ft_meg_source'],'source','-v7.3')
%     save([ssSubjPath(ss) '/ft_meg_freq'],'freq','-v7.3')

end
