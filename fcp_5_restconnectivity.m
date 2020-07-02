function fcp_5_restconnectivity(paths, surf_path, atlas_type)

% FCP_5_RESTCONNECTIVITY estimates functional connectivity using phase-
% locking metrics (PLV, PLI) *or amplitude coupling* (to be added). 
% 
% NOTES:
%   - Will use the fcp_4 subject CSV to pull included participants.  
%
% INPUTS:
%   paths               =   struct defining paths to data, participant
%                           folders, analysis folders, config files, etc. 
%
% OUTPUTS:
%   adjmat:             participant adj matrix preserving results of each
%                       trial, saved into individual analysis folder
%   all_adjmat:         master matrix with all participants adj mats
%                       nodes x nodes x num_participants x num_freq_bands,
%                       saved into group analysis folder
%
% See also: 
%
% Last updated by: Sonya Bells, 2020-02-18
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.
%

%% SETUP

% load config JSON with analysis parameters
config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp4'; %temp for now

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
%         fcp4_output.fig_headsurf  = 'headsurf.png';

%% PARTICIPANT MODELS

rangeOFsubj = 1:length(subj_match.ds);

for ss = rangeOFsubj
     %%% FOR EACH PARTICIPANT --------------------------------------------------
    fprintf('\nWorking on subject %s! \n', subj_match.pid{ss});
 
    %%% LOAD DATA ---------------------------------------------------------
%     load(fullfile(ssSubjPath(ss),sprintf('/%s_meg_sourcemodel_8k_%s',subj_match.pid{ss},'alpha.mat')),'-mat','source');
    load(fullfile(ssSubjPath(ss),sprintf('/%s_meg_sourcemodel_4k_%s',subj_match.pid{ss},'gamma1.mat')),'-mat','source');
   load(fullfile(ssSubjPath(ss),sprintf('/%s_meg_freq_%s',subj_match.pid{ss},'gamma1.mat')),'-mat' ,'freq');
%% compute connectivity
% Compute connectivity matrix between all pairs of dipoles 
% 1) compute imaginary part of the coherency, cfg.method = 'coh'; cfg.complex = 'absimag'
% will return only the imaginary part of the coherence spectrum and effectivly suppress
% spurious coherence driven by electromagnetic field spread (Nolte et al. (2004). Clinical Neurophysiology. 115. 2292-2307)
% source -> contains single trial estimates of amplitude and phase at the source-level (pcc and fourier)
cfg         = [];
cfg.method  ='coh';
cfg.complex = 'absimag';
source_conn = ft_connectivityanalysis(cfg, source);

% figure; imagesc(source_conn.cohspctrm);

%%% Parcellate  -------------------------------
switch atlas_type
    case 'mmp1'
    % using multi-modal parcellation of human cerebral cortex (Glasser et al (2016). Nature. 536. 171)
    %load atlas_MMP1.0_4k.mat - not user specific 
        
        apath      = which('find_megne2.m');
        apath      = apath(1:end-14);
        template    = fullfile(apath, '/functions/functions/atlas_MMP1.0_4k.mat');
        load(template,'-mat', 'atlas');
        atlas.pos = source_conn.pos;

    case 'a2009s'
    %read participants atlas
        mripath = fullfile(surf_path,subj_match.pid{ss});
        datapath = fullfile(mripath,'freesurfer/');
        file_annotation = fullfile(datapath,'/label/lh.aparc.a2009s.annot');
        file_mesh = fullfile(datapath,'/surf/lh.pial');
        atlas = ft_read_atlas({file_annotation, file_mesh},'format','freesurfer_a2009s');

    otherwise 
        ft_error('unsupported format "%s"', atlas_type);

end

% cfg = [];
% cfg.parcellation = 'parcellation';
% cfg.parameter    = 'cohspctrm';
% parc_conn = ft_sourceparcellate(cfg, source_conn, atlas);

cfg = [];
cfg.parcellation = 'parcellation';
cfg.parameter    = 'cohspctrm';
parc_conn = ft_sourceparcellate(cfg, source_conn, atlas);

figure;imagesc(parc_conn.cohspctrm);

%%% Network analysis  -------------------------------
% It is not clear what the effect of the residual spatial leakage of activity is on the estimates of some of these measures
% , so caution should be used interpreting graph metrics derived from connectivity matrices, particularly when 
% comparing groups of experimental participants or experimental conditions
% cfg.method = 'degrees' -> with threshold cfg.threshold = 0.1 -> results in an estimate of the 'node degree' 
% amound of nodes with which a particular node has an estimated connectivity of 0.1 or higher
% 
% cfg           = [];
% cfg.method    = 'degrees';
% cfg.parameter = 'cohspctrm';
% cfg.threshold = .1;
% % network_full = ft_networkanalysis(cfg,source_conn);
% network_parc = ft_networkanalysis(cfg,parc_conn);

% save results
save(fullfile(ssSubjPath(ss),sprintf('/%s_%s_parc',subj_match.pid{ss},'coh_gamma1_a2009s')), 'parc_conn');
% save(fullfile(ssSubjPath(ss),sprintf('/%s_%s_%_parc',subj_match.pid{ss},'coh','degree')), 'network_parc');

%% visualize
% cfg               = [];
% cfg.method        = 'surface';
% cfg.funparameter  = 'degrees';
% cfg.funcolormap   = 'jet';
% ft_sourceplot(cfg, network_full);
% view([-150 30]);
% ft_sourceplot(cfg, network_parc);
% view([-150 30]);

end

