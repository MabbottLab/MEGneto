function fcp_4_Beamforming_Freq_Surf(paths, surf_path)

% fcp_4_Beamforming_Freq_Surf does the following:
% 1) registers CTF mri files with already processed
% freesurfer surf structures and hcp-workbench processing
% 2) frequency analysis on cleaned data
% 3) then carries out beamforming and source projection 
%
%
% NOTES:
%   - Ensure that subj_fcp4.csv is populated with the subject IDs of
%   participants you want to include after checking over initial results. 
%   - Assumes that anatomical MRI is already aligned with MEG data. If not
%   aligned, use ft_volumerealign function. 
%
% INPUTS:
%   paths               =   struct defining paths to data, participant
%                           folders, analysis folders, config files, etc. 
%
% OUTPUTS:
%   .mat
%       ?? - catmatrix:    num_samples x num_trials x num_nodes
%       ?? - srate:        sampling rate
%       ?? - coords:       coordinates of source points w/in
%
% See also: 
%
% Created by: Sonya Bells, 2020-02-14
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.

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

% initialize output files
    % images
        fcp1_output.fig_headsurf  = 'headsurf.png';

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

%     %% VISUALIZATION: Headmodel and Surfaces --------------------------------
%     figure;
%         % make the headmodel surface transparent
%         ft_plot_headmodel(headmodel, 'edgecolor', 'none'); alpha 0.9
%         ft_plot_mesh(ft_convert_units(sourcemodel_4k, 'cm'),'vertexcolor',sourcemodel_4k.sulc);
%         ft_plot_sens(data.grad);
%         view([0 -90 0])

    %% Prepare for beamforming
    
    sourcemodel_4k = ft_convert_units(sourcemodel_4k, 'cm');
    %%% Compute the leadfield ----------------------------------------------------
    cfg             = [];
    cfg.grid        = sourcemodel_4k;
    cfg.headmodel   = headmodel;
    % cfg.channel     = {'MEG'};
    lf              = ft_prepare_leadfield(cfg, data);

    %%% Spectral analysis ---------------------------------------------------------
    % Compute sensor level Fourier spectra, to be used for cross-spectral density computation
    cfg             = [];
    cfg.method      = config.beamforming.freqDomain.method; %'mtmfft'; % analysis an entire spectrum for entire data length (mulitaper freq transformation)
    cfg.output      = config.beamforming.freqDomain.output; %'fourier'; % return complex Fourier-spcetra ; powandcsd return the power and cross-spectra
    cfg.taper       = config.beamforming.freqDomain.taper; %'dpss'; %default  
    cfg.tapsmofrq   = config.beamforming.freqDomain.tapsmofrq; %1;
    cfg.foi         = config.beamforming.freqDomain.foi; %10
    cfg.keeptrials  = config.beamforming.options.keeptrials; %'yes'
    
    % check whether conn method is wPLI
    if contains(config.connectivity.method, "wPLI")
        cfg.t_ftimwin = ones(length(cfg.foi),1).*0.5; %length of time window in seconds = 0.5s
        cfg.toi       = 0:0.05:1.5; %the times of which the analysis windows should be centered (in seconds) from 0 to 6.99s in steps of 0.05
    end
    
    % actual frequency analysis
    freq           = ft_freqanalysis(cfg, data);

    %%% Forward Model ---------------------------------------------------------
    % Use a cortical sheeet based source model, in which individual dipole locations 
    % are constrated to cortical sheet. 

    % Source reconstructrion    
    cfg                 = [];
    cfg.frequency       = freq.freq;
    cfg.method          = config.beamforming.method; % eloreta or pcc (coh)
    cfg.grid            = lf;
    cfg.headmodel       = headmodel;
    cfg.keeptrials      = 'yes';
    cfg.lambda          = '10%';
    cfg.projectnoise    = 'yes';
    cfg.fixedori        = 'yes';  
    source              = ft_sourceanalysis(cfg, freq);

    % save results - Need to modify
    save(fullfile(ssSubjPath(ss),sprintf('/%s_meg_sourcemodel_4k_%s',subj_match.pid{ss},'gamma1')), 'source');
    save(fullfile(ssSubjPath(ss),sprintf('/%s_meg_freq_%s',subj_match.pid{ss},'gamma1')), 'freq');

    %save([ssSubjPath(ss) '/ft_meg_source'],'source','-v7.3')
    %save([ssSubjPath(ss) '/ft_meg_freq'],'freq','-v7.3')

end
