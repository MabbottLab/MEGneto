function fcp_5_restconnectivity(paths)

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

%figure; imagesc(source_conn.cohspctrm);

%%% Parcellate  -------------------------------
% using multi-modal parcellation of human cerebral cortex (Glasser et al (2016). Nature. 536. 171)
load atlas_MMP1.0_4k.mat; 
atlas.pos = source_conn.pos;

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

cfg           = [];
cfg.method    = 'degrees';
cfg.parameter = 'cohspctrm';
cfg.threshold = .1;
network_full = ft_networkanalysis(cfg,source_conn);
network_parc = ft_networkanalysis(cfg,parc_conn);
%% visualize
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'degrees';
cfg.funcolormap   = 'jet';
ft_sourceplot(cfg, network_full);
view([-150 30]);
ft_sourceplot(cfg, network_parc);
view([-150 30]);

