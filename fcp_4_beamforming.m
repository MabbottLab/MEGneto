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
template_seg = ft_volumesegment(cfg, template); % unit: mm

%%% PREPARE HEAD MODEL WITH SEGMENTED TEMPLATE BRAIN ----------------------
cfg                 = [];
cfg.method          = config.beamforming.headmodel.method; 
template_headmodel  = ft_prepare_headmodel(cfg, template_seg); % unit: mm
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
%   figure
%   hold on
%   ft_plot_vol(template_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
%   ft_plot_mesh(template_grid.pos(template_grid.inside,:));
% 
% % VISUALIZATION: ALIGNMENT ----------------------------------------------

% %% load atlas and create a binary mask
%   atlas = ft_read_atlas([ftpath '/' config.beamforming.atlas.filepath]);
%   atlas = ft_convert_units(atlas, config.beamforming.headmodel.units);% assure that atlas and template_grid are expressed in the %same units
% 
% %% get internal tissue labels
%   cfg         = [];
%   cfg.atlas   = atlas;
%   cfg.roi     = atlas.tissuelabel;
%   cfg.inputcoord = config.beamforming.atlas.inputcoord;
%   mask        = ft_volumelookup(cfg,template_grid);
% 
% %% create temporary mask according to the atlas entries
%   tmp           = repmat(template_grid.inside,1,1);
%   tmp(tmp==1)   = 0;
%   tmp(mask)     = 1;
% 
% %% define inside locations according to the atlas based mask
%   template_grid.inside = tmp;
% 
% %% plot the atlas based grid
%   figure;
%   ft_plot_mesh(template_grid.pos(template_grid.inside,:));

%% PARTICIPANT MODELS

rangeOFsubj = 1:length(subj_match.ds);

for ss = rangeOFsubj
%%% FOR EACH PARTICIPANT --------------------------------------------------
    fprintf('\nWorking on subject %s! \n', subj_match.pid{ss});

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
%     cfg.title           = ['Segmentation: ', subj_match.pid{ss}];
%     ft_sourceplot(cfg, seg);
%  
%     hf = gcf;
%     hf.Position(1:2) = [10 10]; 
%     hf.Position(3:4) = (800 / hf.Position(4)) .* hf.Position(3:4);
%     print(hf, [ssSubjPath(ss) '/segmented_mri_alignment'], '-dpng', '-r600');
    
%%% PREPARE SUBJECT-SPECIFIC HEAD MODEL WITH TEMPLATE HEAD MODEL ----------
    cfg                = [];
    cfg.grid.warpmni   = config.beamforming.subj.grid.warpmni;
    cfg.grid.template  = template_grid;
    cfg.grid.nonlinear = config.beamforming.subj.grid.nonlinear;
    cfg.mri            = mri;
    cfg.grid.unit      = config.beamforming.subj.grid.unit; %mm
    grid               = ft_prepare_sourcemodel(cfg);
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
%     alpha 0.4; % opacity of headmodel
%     ft_plot_mesh(grid.pos(grid.inside,:)); % plot source model
%     ft_plot_sens(data.grad,'style','ob'); % plot MEG channels (sensors)
%     hold off;
%     view(45,10);
% 
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
    
    %%% project all trials through common spatial filter
    cfg                 = [];
    cfg.grid             = grid; % source model
    cfg.grid.filter      = source_t_avg.avg.filter;
    cfg.grid.leadfield   = leadfield.leadfield;
    cfg.grad             = data.grad; % sensor position (gradiometer)
    cfg.headmodel        = hdm;
    cfg.method           = config.beamforming.method;
    cfg.keeptrials       = 'yes';
    % cfg.keepfilter       = config.beamforming.options.keepfilter;
    cfg.rawtrial         = config.beamforming.options.rawtrial;
    source_t_trials      = ft_sourceanalysis(cfg, tlock);
    
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
    fullPath        = which('ft_preprocessing.m');
    [pathstr,~,~]   = fileparts(fullPath);
    atlas           = ft_read_atlas([pathstr, '/template/atlas/aal/ROI_MNI_V4.nii']);
    atlas           = ft_convert_units(atlas, 'cm');

    % source interpolate
    cfg              = [];
    cfg.interpmethod = 'nearest';
    cfg.parameter    = 'tissue';
    source_aal       = ft_sourceinterpolate(cfg,atlas,sourcemodel);

    % actual interpolation
    vector_virtual_sources      = [];   % overall AAL region timeseries across trial
    ori_avg                     = [];   % orientation average

    %%% FOR EACH TRIAL ----------------------------------------------------
    for t = 1:projection.df
        vector_virtual_sources_trial = [];
        ori_avg_trial                = [];
        ori_116                      = [];
        %%% AND FOR EACH AAL REGION ---------------------------------------
        for i = 1:90
            source_timeseries        = [];
            ori_region               = [];
            
            % identify source coords that fall within AAL region
            aal_node                 = find(source_aal.tissue==i); 
            sources_innode           = projection.trial(t).mom(aal_node); 
            ori_innode               = projection.trial(t).ori(aal_node); 

            % isolate source timeseries in region
            for j = 1:length(sources_innode)    %%% FOR EACH SOURCE W/IN THE AAL REGION
                source_timeseries   = cat(1, source_timeseries, sources_innode{j}); 
                ori_region          = cat(1, ori_region, ori_innode{1,j});
            end
            
            % IF MORE THAN 1 SOURCE POINT W/IN NODE
            if length(sources_innode) > 1
                vector_virtual_sources_trial = ...
                    cat(1, vector_virtual_sources_trial, mean(source_timeseries));
                ori_avg_trial                = ...
                    cat(1,ori_avg_trial,mean(ori_region));      
            % IF NO SOURCE POINTS W/IN NODE
            elseif isempty(sources_innode)
                warning('NO NODE %d\n',i);
            % IF ONLY 1 SOURCE POINT W/IN NODE
            else
                vector_virtual_sources_trial = ...
                    cat(1,vector_virtual_sources_trial,source_timeseries);
                ori_avg_trial                = ...
                    cat(1,ori_avg_trial,ori_region);
            end
            
            % orientation of eigenvectors of each source (divided into AAL regions)
            ori_region_padded = zeros(ceil(0.03*nnz(source_aal.tissue)),3); 
            ori_region_padded(1:size(ori_region,1),:) = ori_region;
            ori_116 = cat(3,ori_116,ori_region_padded); 
        end
        
        %%% CONCATENATE ACROSS TRIALS -------------------------------------
        vector_virtual_sources = cat(3,vector_virtual_sources,vector_virtual_sources_trial);
        ori_avg                = cat(3,ori_avg,ori_avg_trial); 
    end
    
%%% FINAL PARTICIPANT OUTPUT VARIABLES ------------------------------------
    catmatrix = permute(vector_virtual_sources,[2 3 1]); % reorder dimensions for next processing steps
    srate     = data.fsample;                            % keep track of sampling rate
    coords    = projection.pos;                          % coordinates

%%% SAVE OUTPUT -----------------------------------------------------------
    save([ssSubjPath(ss) '/AAL_beamforming_results'],'catmatrix','srate','coords','-mat','-v7.3')

%%% OPTIMIZING RUN SPACE --------------------------------------------------
    clear coords catmatrix srate vector_virtual_sources vector_virtual_sources_trial source_timeseries ...
        sources_innode ori_avg ori_avg_trial ori_116 atlas sourcemodel source_t_trials projection seg mri data grid hdm ...
        source_t_avg tlock leadfield aal_node ori_region ori_innode ori_region_padded
    fprintf('\nDone subject %s! \n', subj_match.pid{ss});
end
