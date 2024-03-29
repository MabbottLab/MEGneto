function fcp_4_beamforming_tmplGrid(paths)

% FCP_4_BEAMFORMING carries out beamforming and source projection on
% cleaned data. A T1 template head model is loaded in (normalizing across
% particpants) to create a volume conduction model and subsequently, a
% source model. Each participant's MRI data is then loaded in and their
% head and source models are created. After a leadfield is computed for
% efficienct inverse modelling, beamforming is performed by computing a
% covariance matrix and projecting each participant's trials through a
% spatial filter. Once the virtual sources are projected to their dominant
% orientation, the source data is interpolated onto an atlas. All
% reconstructed sources that fall into regions of interest on the atlas,
% are averaged to generate a representative timseries for the region.
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
%       4. MMP atlas (Glasser, 2016)
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
% See also: DS_PID_MATCH, WRITE_MATCH_IF_NOT_EMPTY, FT_READ_MRI,
% FT_VOLUMESEGMENT, FT_CONVERT_UNITS, FT_PREPARE_HEADMODEL,
% FT_PREPARE_SOURCEMODEL, FT_RESAMPLEDATA, FT_PREPARE_LEADFIELD,
% FT_TIMELOCKANALYSIS, FT_SOURCEANALYSIS, FT_SOURCEDESCRIPTIVES,
% FT_READ_ATLAS, FT_SOURCEINTERPOLATE, FT_VOLUMELOOKUP
%
% Last updated by: Julie Tseng, 2020-01-08
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.

%% SET UP LOGGING FILE

right_now = clock;
log_filename = [paths.conf_dir '/log_' sprintf('%02.f%02.f%02.f', right_now(1:3))];
diary(log_filename)

fprintf('\n\n%02.f:%02.f:%02.f       Now running **%s**.\n', ...
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

ft_path = fileparts(which('ft_defaults.m'));
template_grid       = load([ft_path '/template/sourcemodel/standard_sourcemodel3d7point5mm.mat']);
template_grid = template_grid.sourcemodel;

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

for ss = rangeOFsubj % for each participant that has matched MEG/MRI data
%%% FOR EACH PARTICIPANT --------------------------------------------------
    right_now = clock;
    fprintf('%02.f:%02.f:%02.f      Working on subject %s!\n', ...
        right_now(4:6), subj_match.pid{ss})

%%% LOAD ANATOMICAL MRI DATA ----------------------------------------------
    mri     = ft_read_mri([paths.rawmri '/' subj_match.pid{ss} '_V2.mri']);
    mri     = ft_convert_units(mri,'cm');
    
    % check for fiducials which help to localize head position relative to
    % the sensors
    if any(mri.hdr.fiducial.mri.nas) == 0 || any(mri.hdr.fiducial.mri.lpa) == 0  || any(mri.hdr.fiducial.mri.rpa) == 0
        error('No fiducials found for subject %s!', subj_match.pid{ss});
    end

%%% SEGMENT ANATOMICAL MRI ------------------------------------------------
    cfg        = []; % set up config for volume segmentation
    cfg.output = 'brain';
    seg        = ft_volumesegment(cfg, mri); % segment participant MRI

%%% PREPARE HEAD MODEL WITH SEGMENTED PARTICIPANT BRAIN -------------------
    cfg             = []; % set up confirm to prepare the participant head model
    cfg.method      = config.beamforming.headmodel.method;
    headmodel_units = 'cm';
    hdm             = ft_prepare_headmodel(cfg, seg); % prepare head model
    hdm             = ft_convert_units(hdm,headmodel_units); % convert to specified units
    
%%% LOAD MEG DATA ---------------------------------------------------------
    load([ssSubjPath(ss) '/ft_meg_fullyProcessed.mat'],'-mat','data');
    
    % Resample all datasets
    if ~(config.filteringParameters.sampleRate == data.fsample)
        cfg             = []; % set up config to resample the data
        cfg.resamplefs  = config.filteringParameters.sampleRate;
        cfg.detrend     = 'no';
        data_resamp     = ft_resampledata(cfg, data); % resample data to specified rate
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
    
%%% PREPARE SUBJECT-SPECIFIC SOURCE MODEL WITH TEMPLATE HEAD MODEL ----------
    cfg                 = []; % set upp config for participant source model preparation 
    cfg.template        = template_grid;
    cfg.warpmni         = config.beamforming.subj.grid.warpmni;
    cfg.nonlinear       = config.beamforming.subj.grid.nonlinear;
    cfg.mri             = mri;
    cfg.unit            = config.beamforming.subj.grid.unit; 
    grid                = ft_prepare_sourcemodel(cfg); % prepare source model
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
    % the leadfield is used to provide information on the contribution of a
    % dipole source at a given location in a sensor's region
    cfg              = []; % set up config to prepare the leadfield
    cfg.headmodel    = hdm;
    cfg.channel      = 'MEG';
    cfg.reducerank   = 2;
    cfg.grid.pos     = grid.pos;
    cfg.grid.inside  = grid.inside;
    cfg.normalize    = config.beamforming.leadfield.normalize;
    cfg.grad         = data.grad;
    leadfield        = ft_prepare_leadfield(cfg, data); % create leadfield

%% ACTUAL BEAMFORMING

    right_now = clock;
    fprintf('%02.f:%02.f:%02.f       Now, for actual beamforming...\n', ...
        right_now(4:6))

%%% VECTOR - Time Domain Source Reconstruction ----------------------------

    %%% compute common spatial filter (returns: COVARIANCE MATRIX)
    % the covariance matrix tells us how related the sensors are
    cfg                    = []; % set up config to compute covariance matrix
    cfg.covariance         = config.beamforming.timeDomain.covariance;
    cfg.channel            = 'MEG';
    cfg.covariancewindow   = config.beamforming.timeDomain.covariancewindow;
    cfg.vartrllength       = config.beamforming.timeDomain.vartrllength;
    cfg.keeptrials         = config.beamforming.options.keeptrials;
    tlock                  = ft_timelockanalysis(cfg, data); % compute covariance matrix

    %%% calculate sensor weights (actual beamforming)
    cfg                 = []; % set up config for sensor weight calculation
    cfg.grad            = data.grad;            % sensor position (gradiometer)
    cfg.headmodel       = hdm;
    cfg.grid            = grid;          % source model
    cfg.grid.leadfield  = leadfield.leadfield;
    cfg.keepfilter      = config.beamforming.options.keepfilter;
    source_t_avg        = ft_sourceanalysis(cfg, tlock);
    
    %%% project all trials thru spatial filter
    cfg                 = []; % set up config for beamforming
    cfg.grid             = grid; % source model
    cfg.grid.filter      = source_t_avg.avg.filter;
    cfg.grid.leadfield   = leadfield.leadfield;
    cfg.grad             = data.grad; % sensor position (gradiometer)
    cfg.headmodel        = hdm;
    cfg.method           = config.beamforming.method;
    cfg.keeptrials       = 'yes';
    cfg.rawtrial         = 'yes';       
    source_t_trials      = ft_sourceanalysis(cfg, tlock); % perform beamforming
        
    %%% project virtual sources to strongest (dominant) orientation
    %%% (taking the largest eigenvector of the sources timeseries)
    cfg                  = []; % set up config for projecting to dominant orientation
    cfg.projectmom       = config.beamforming.timeDomain.projectmom;
    cfg.keeptrials       = 'yes';
    projection           = ft_sourcedescriptives(cfg, source_t_trials); % project to dominant orientation
    
%%% INTERPOLATE ATLAS ONTO VIRTUAL SOURCES ----------------------------

    % setup for atlas interpolation - get coordinates
    sourcemodel.pos = template_grid.pos; 
    
    % load atlas
    if contains(config.beamforming.atlas.filepath, 'mmp') % if MMP glasser atlas
        megneto_path        = fileparts(which('fcp_4_beamforming.m'));
        atlas               = ft_read_atlas([megneto_path '/external/atlas/mmp.mat']);
    else
        fullPath                = which('ft_preprocessing.m');
        [pathstr,~,~]           = fileparts(fullPath);
        atlas                   = ft_read_atlas([pathstr config.beamforming.atlas.filepath]);
        if contains(config.beamforming.atlas.filepath, 'aal')
            atlas.tissuelabel   = atlas.tissuelabel(1:90); % we only want non-cerebellar regions (isolate desired regions)
            atlas.tissue(atlas.tissue > 90) = 0;
        end
    end
    atlas           = ft_convert_units(atlas, 'cm'); % convert atlas units to centimeters

    % source interpolate
    cfg              = [];
    cfg.interpmethod = 'nearest';
    cfg.parameter    = 'tissue';
    source_atlas       = ft_sourceinterpolate(cfg,atlas,sourcemodel); % interpolate source activity onto voxels of anatomical description of the brain

    % actual interpolation
    catmatrix      = NaN(length(projection.time), ... % set up empty matrix to store reconstructed timeseries for each trial and region of interest
                         length(projection.trial), ...
                         length(source_atlas.tissuelabel));   % over all ROI timeseries across trial

    var_explained  = NaN(1, ... % set up empty matrix to store variance explained of first principial component for each trial and region of interest 
                         length(projection.trial), ...
                         length(source_atlas.tissuelabel));
    %%% FOR EACH TRIAL ----------------------------------------------------
    right_now = clock;
    fprintf('%02.f:%02.f:%02.f       Identifying ROI timeseries!\n', ...
        right_now(4:6))
    
    for t = 1:projection.df
        %%% AND FOR EACH NODE ---------------------------------------------
        for i = 1:max(size(atlas.tissuelabel))
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
                % ori_avg(:,t,i) = nanmean(ori_region,1);
            % IF NO SOURCE POINTS W/IN NODE
            else
                warning('No sources in ROI %s.\n',atlas.tissuelabel{i});
            end
        end
    end
    
%%% FINAL PARTICIPANT OUTPUT VARIABLES ------------------------------------
    srate     = data.fsample;                            % keep track of sampling rate
    coords    = projection.pos;                          % coordinates

%%% SAVE OUTPUT -----------------------------------------------------------
    if isnan(var_explained) % if the variance explained has not been populated, don't save it
        save([ssSubjPath(ss) '/atlas_beamforming_results.mat'],'catmatrix', 'srate','coords','-mat','-v7.3')
    else
        save(['/home/kwalonso/TP/' subj_match.pid{ss} '/atlas_beamforming_results.mat'],'catmatrix', 'var_explained', 'srate','coords','-mat','-v7.3')
    end 
        
%%% OPTIMIZING RUN SPACE --------------------------------------------------
    clear coords catmatrix srate source_timeseries ...
        atlas sourcemodel source_t_trials projection seg mri data grid hdm ...
        source_t_avg tlock leadfield 
    
    right_now = clock;
    fprintf('%02.f:%02.f:%02.f       Done subject %s!\n', ...
        right_now(4:6), subj_match.pid{ss})
    
end

%% turn off diary
right_now = clock;
fprintf('%02.f:%02.f:%02.f       Done running **%s**.\n', ...
    right_now(4:6), mfilename)

%% send email to user
sendEmail("beamforming", string(config.contact));

diary off
