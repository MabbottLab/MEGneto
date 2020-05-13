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
%   - Available atlases:
%       1. AAL 116 atlas (all areas <=90 to exclude cerebellum)
%       2. Brainnetome atlas
%       3. Yeo atlas (7 network or 17 network)
%
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

%% TEMPLATE HEAD MODEL/SOURCE MODEL

%%% LOAD T1 TEMPLATE FROM SPM W/IN FIELDTRIP ------------------------------
ftpath      = which('ft_defaults.m');
ftpath      = ftpath(1:end-14);
template    = ft_read_mri(fullfile(ftpath, '/external/spm8/templates/T1.nii'));
template.coordsys = config.beamforming.template.coordsys;

%%% SEGMENT TEMPLATE BRAIN ------------------------------------------------
    % needed to describe boundaries that define which dipole locations are 
    % 'inside' the brain.
cfg          = [];
template_seg = ft_volumesegment(cfg, template); 
template_seg = ft_convert_units(template_seg, 'cm');

%%% PREPARE HEAD MODEL WITH SEGMENTED TEMPLATE BRAIN ----------------------
cfg                 = [];
cfg.method          = config.beamforming.headmodel.method; 
template_headmodel  = ft_prepare_headmodel(cfg, template_seg); 
template_headmodel  = ft_convert_units(template_headmodel, ...
                        config.beamforming.headmodel.units); 
%%% CONSTRUCT DIPOLE GRID IN TEMPLATE BRAIN COORDINATES -------------------               
cfg                 = [];
cfg.resolution      = config.beamforming.template.grid.resolution;
cfg.tight           = config.beamforming.template.grid.tight;
cfg.inwardshift     = config.beamforming.template.grid.inwardshift; % negative for inside
cfg.headmodel       = template_headmodel;
template_grid       = ft_prepare_sourcemodel(cfg);

%%% VISUALIZATION: DISPLAY HEAD MODEL -------------------------------------
%       figure
%       hold on
%       ft_plot_vol(template_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
%       ft_plot_mesh(template_grid.pos(template_grid.inside,:));
% 
% % VISUALIZATION: ALIGNMENT ----------------------------------------------

% %% load atlas and create a binary mask
%   atlas = ft_read_atlas([ftpath '/' config.beamforming.atlas.filepath]);
%   atlas = ft_convert_units(atlas, config.beamforming.headmodel.units);% assure that atlas and template_grid are expressed in the %same units
% % 
% %%% get internal tissue labels
%   cfg         = [];
%   cfg.atlas   = atlas;
%   cfg.roi     = atlas.tissuelabel;
%   cfg.inputcoord = config.beamforming.atlas.inputcoord;
%   mask        = ft_volumelookup(cfg,template_grid);
% % 
% % %% create temporary mask according to the atlas entries
%   tmp           = repmat(template_grid.inside,1,1);
%   tmp(tmp==1)   = 0;
%   tmp(mask)     = 1;
% 
% %%% define inside locations according to the atlas based mask
%   template_grid.inside = tmp;
% 
% %%% plot the atlas based grid
%   figure;
%   ft_plot_mesh(template_grid.pos(template_grid.inside,:));

%% PARTICIPANT MODELS

rangeOFsubj = 1:length(subj_match.ds);

for ss = rangeOFsubj
%%% FOR EACH PARTICIPANT --------------------------------------------------
    right_now = clock;
    fprintf('%d:%d:%02.f       Working on subject %s!\n', ...
        right_now(4:6), subj_match.pid{ss})

%%% LOAD ANATOMICAL MRI DATA ----------------------------------------------
    mri     = ft_read_mri([paths.rawmri '/' subj_match.pid{ss} '_V2.mri']);
    mri     = ft_convert_units(mri,'cm');
    
    % check for fiducials
    if any(mri.hdr.fiducial.mri.nas) == 0 || any(mri.hdr.fiducial.mri.lpa) == 0  || any(mri.hdr.fiducial.mri.rpa) == 0
        error('No fiducials found for subject %s!', subj_match.pid{ss});
    end

%%% SEGMENT ANATOMICAL MRI ------------------------------------------------
    cfg        = [];
    cfg.output = 'brain';
    seg        = ft_volumesegment(cfg, mri);

%%% PREPARE HEAD MODEL WITH SEGMENTED PARTICIPANT BRAIN -------------------
    cfg             = [];
    cfg.method      = config.beamforming.headmodel.method;
    headmodel_units = 'cm';
    hdm             = ft_prepare_headmodel(cfg, seg);
    hdm             = ft_convert_units(hdm,headmodel_units);
    
%%% LOAD MEG DATA ---------------------------------------------------------
    load([ssSubjPath(ss) '/ft_meg_fullyProcessed.mat'],'-mat','data');
    
    % Resample all datasets
    if ~(config.filteringParameters.sampleRate == data.fsample)
        cfg             = [];
        cfg.resamplefs  = config.filteringParameters.sampleRate;
        cfg.detrend     = 'no';
        data_resamp     = ft_resampledata(cfg, data);
        data = data_resamp;
        clear data_resamp
    end

%%% VISUALIZATION: CHECK SEGMENTED MRI FOR PROPER ALIGNMENT ---------------
%     seg.transform  = mri.transform;
%     seg.anatomy    = mri.anatomy;
% 
%     cfg = [];
%     cfg.method          = config.beamforming.checkMRIvolumes.method;
%     cfg.slicesdim       = config.beamforming.checkMRIvolumes.slicesdim;
%     cfg.nslices         = config.beamforming.checkMRIvolumes.nslices;
%     cfg.anaparameter    = 0.5;
%     cfg.funparameter    = 'brain';
%     cfg.title           = ['Segmentation: ', subj_match.pid{ss}];
%     ft_sourceplot(cfg, seg);
%  
%     hf = gcf;
%     hf.Position(1:2) = [10 10]; 
%     hf.Position(3:4) = (800 / hf.Position(4)) .* hf.Position(3:4);
%     print(hf, [ssSubjPath(ss) '/segmented_mri_alignment'], '-dpng', '-r600');
    
%%% PREPARE SUBJECT-SPECIFIC HEAD MODEL WITH TEMPLATE HEAD MODEL ----------
    cfg                 = [];
    cfg.template        = template_grid;
    cfg.warpmni         = config.beamforming.subj.grid.warpmni;
    cfg.nonlinear       = config.beamforming.subj.grid.nonlinear;
    cfg.mri             = mri;
    cfg.unit            = config.beamforming.subj.grid.unit; 
    grid                = ft_prepare_sourcemodel(cfg);
%    coords             = sourcemodel.pos;

%%% VISUALIZATION: TISSUE PROBABILITY MAPS --------------------------------
%%% make a figure of the single subject headmodel, and grid positions
%     figure; hold on;
%     ft_plot_vol(hdm, 'edgecolor', 'none', 'facealpha', 0.4);
%     ft_plot_mesh(grid.pos(grid.inside,:));

%%% check alignment of source model and head model and save as image
%     figure
%     hold on;
%     ft_plot_vol(hdm,'edgecolor','none','facecolor', 'cortex'); % plot head model
%     alpha 0.9; % opacity of headmodel
%     ft_plot_mesh(grid.pos(grid.inside,:)); % plot source model
%     ft_plot_sens(data.grad,'style','ob'); % plot MEG channels (sensors)
%     hold off;
%     view(45,10);

%     hf = gcf;
%     hf.Position(1:2) = [10 10];
%     hf.Position(3:4) = (800 / hf.Position(4)) .* hf.Position(3:4);
%     print(hf, [ssSubjPath(ss) '/source_head_align'], '-dpng', '-r600');

%%% COMPUTE LEADFIELD -----------------------------------------------------
    cfg              = [];
    cfg.headmodel    = hdm;
    cfg.channel      = 'MEG';
    cfg.reducerank   = 2;
    cfg.grid.pos     = grid.pos;
    cfg.grid.inside  = grid.inside;
    cfg.normalize    = config.beamforming.leadfield.normalize;
    cfg.grad         = data.grad;
    leadfield        = ft_prepare_leadfield(cfg, data);

%% ACTUAL BEAMFORMING

    right_now = clock;
    fprintf('%d:%d:%02.f       Now, for actual beamforming...\n', ...
        right_now(4:6))

%%% VECTOR - Time Domain Source Reconstruction ----------------------------

    %%% compute common spatial filter (returns: COVARIANCE MATRIX)
    cfg                    = [];
    cfg.covariance         = config.beamforming.timeDomain.covariance;
    cfg.channel            = 'MEG';
    cfg.covariancewindow   = config.beamforming.timeDomain.covariancewindow;
    cfg.vartrllength       = config.beamforming.timeDomain.vartrllength;
    cfg.keeptrials         = config.beamforming.options.keeptrials;
    tlock                  = ft_timelockanalysis(cfg, data);  

    %%% calculate sensor weights (actual beamforming)
    cfg                 = [];
    cfg.grad            = data.grad;            % sensor position (gradiometer)
    cfg.headmodel       = hdm;
    cfg.grid            = grid;          % source model
    cfg.grid.leadfield  = leadfield.leadfield;
    cfg.keepfilter      = config.beamforming.options.keepfilter;
    source_t_avg        = ft_sourceanalysis(cfg, tlock);
    
    start_idx = 1:100:length(tlock.trialinfo);
    end_idx = start_idx + 99;
    end_idx(end) = length(tlock.trialinfo);
    
    for part = 1:(ceil(length(tlock.trialinfo)/100))
        %%% project all trials through common spatial filter
        cfg                 = [];
        cfg.grid             = grid; % source model
        cfg.grid.filter      = source_t_avg.avg.filter;
        cfg.grid.leadfield   = leadfield.leadfield;
        cfg.grad             = data.grad; % sensor position (gradiometer)
        cfg.headmodel        = hdm;
        cfg.method           = config.beamforming.method;
        cfg.keeptrials       = 'yes';
        cfg.rawtrial         = config.beamforming.options.rawtrial;
        
        cfg2                 = [];
        cfg2.trials          = (start_idx(part):end_idx(part));
        
        source_t_trials      = ft_sourceanalysis(cfg, ft_selectdata(cfg2, tlock));
        
        %%% project virtual sources to strongest (dominant) orientation
        %%% (taking the largest eigenvector of the sources timeseries)
        cfg                  = [];
        cfg.projectmom       = config.beamforming.timeDomain.projectmom;
        cfg.keeptrials       = 'yes';
        projection{part}     = ft_sourcedescriptives(cfg, source_t_trials);
    end
    
    %%% project virtual sources to strongest (dominant) orientation
    %%% (taking the largest eigenvector of the sources timeseries)
    cfg                  = [];
    cfg.projectmom       = config.beamforming.timeDomain.projectmom;
    cfg.keeptrials       = 'yes';
    projection           = ft_sourcedescriptives(cfg, source_t_trials);
    
%%% INTERPOLATE AAL ATLAS ONTO VIRTUAL SOURCES ----------------------------

    % setup for AAL interpolation - get coordinates
    sourcemodel.pos = template_grid.pos; 
    
    % load atlas
    fullPath                = which('ft_preprocessing.m');
    [pathstr,~,~]           = fileparts(fullPath);
    atlas                   = ft_read_atlas([pathstr config.beamforming.atlas.filepath]);
    if contains(config.beamforming.atlas.filepath, 'aal')
        atlas.tissuelabel   = atlas.tissuelabel(1:90); % we only want non-cerebellar regions
        atlas.tissue(atlas.tissue > 90) = 0;
    end
    atlas           = ft_convert_units(atlas, 'cm');

    % source interpolate
    cfg              = [];
    cfg.interpmethod = 'nearest';
    cfg.parameter    = 'tissue';
    source_atlas       = ft_sourceinterpolate(cfg,atlas,sourcemodel);

    % actual interpolation
    catmatrix      = NaN(length(projection.time), ...
                         length(projection.trial), ...
                         length(source_atlas.tissuelabel));   % overall AAL region timeseries across trial

    %%% FOR EACH TRIAL ----------------------------------------------------
    right_now = clock;
    fprintf('%d:%d:%02.f       Projecting to AAL sources!\n', ...
        right_now(4:6))
    
    for t = 1:projection.df
        %%% AND FOR EACH NODE ---------------------------------------------
        for i = 1:max(size(atlas.tissuelabel))
            % identify source coords that fall within AAL region
            node                     = find(source_atlas.tissue==i); 
            source_timeseries        = cell2mat(projection.trial(t).mom(node)'); % get the timeseries; num_nodes x time
            ori_region               = cell2mat(projection.trial(t).ori(node)'); % orientations; num_nodes x time
            
            % IF NODE EXISTS
            if size(source_timeseries, 1) >= 1
                catmatrix(:,t,i) = nanmean(source_timeseries,1); % take avg across source points
                ori_avg(:,t,i) = nanmean(ori_region,1);
            % IF NO SOURCE POINTS W/IN NODE
            else
                warning('NO NODE %d\n',i);
            end
        end
    end
    
%%% FINAL PARTICIPANT OUTPUT VARIABLES ------------------------------------
    srate     = data.fsample;                            % keep track of sampling rate
    coords    = projection.pos;                          % coordinates

%%% SAVE OUTPUT -----------------------------------------------------------
    save([ssSubjPath(ss) '/AAL_beamforming_results'],'catmatrix','srate','coords','-mat','-v7.3')

%%% OPTIMIZING RUN SPACE --------------------------------------------------
    clear coords catmatrix srate source_timeseries ...
        atlas sourcemodel source_t_trials projection seg mri data grid hdm ...
        source_t_avg tlock leadfield aal_node
    
    right_now = clock;
    fprintf('%d:%d:%02.f       Done subject %s!\n', ...
        right_now(4:6), subj_match.pid{ss})
    
end

%% turn off diary
right_now = clock;
fprintf('%d:%d:%02.f       Done running **%s**.\n', ...
    right_now(4:6), mfilename)

%% send email to user
sendEmail("beamforming", string(config.contact));

diary off
