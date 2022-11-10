function Step3_Beamforming(config, pid)

% to fill this part in...

%% SET UP LOGGING FILE

this_output = [config.meta.project_path '/' config.meta.analysis_name '/' pid];
load([this_output '/out_struct.mat'])

%% TEMPLATE SOURCE MODEL

% load from FieldTrip templates
load(sprintf('%s/template/sourcemodel/standard_sourcemodel3d%dmm', ...
                config.meta.fieldtrip_path, ...
                config.step3.templateRes), 'sourcemodel');
template_grid = ft_convert_units(sourcemodel, 'mm');
clear sourcemodel;

%% PARTICIPANT MODELS

%%% LOAD ANATOMICAL MRI DATA ----------------------------------------------
    load([this_output '/Prep_T1_aligned.mat']);
    mri     = ft_convert_units(mri,'cm');
    
    % check for fiducials which help to localize head position relative to
    % the sensors
    if any(mri.cfg.fiducial.nas) == 0 || any(mri.cfg.fiducial.lpa) == 0  || any(mri.cfg.fiducial.rpa) == 0
        error('No fiducials found for subject %s!', pid);
    end

%%% SEGMENT ANATOMICAL MRI ------------------------------------------------
    cfg        = []; % set up config for volume segmentation
    cfg.output = 'brain';
    seg        = ft_volumesegment(cfg, mri); % segment participant MRI

%%% PREPARE HEAD MODEL WITH SEGMENTED PARTICIPANT BRAIN -------------------
    cfg             = []; % set up confirm to prepare the participant head model
    cfg.method      = 'singleshell';
    hdm             = ft_prepare_headmodel(cfg, seg); % prepare head model
    
%%% LOAD MEG DATA ---------------------------------------------------------
    load([this_output '/step2_data_fullyProcessed.mat']);
    data = data_clean;
    clear data_clean
    
%%% PREPARE SUBJECT-SPECIFIC SOURCE MODEL WITH TEMPLATE HEAD MODEL ----------
    cfg                 = []; % set upp config for participant source model preparation 
    cfg.mri             = mri;
    cfg.nonlinear       = 'yes';
    cfg.unit            = 'mm';
    cfg.template        = template_grid;
    cfg.spmversion      = 'spm8';
    cfg.method          = 'basedonmni';
    grid                = ft_prepare_sourcemodel(cfg); % prepare source model

%%% VISUALIZATION
%%% check alignment of source model and head model and save as image
    figure
    hold on;
    ft_plot_headmodel(hdm,'edgecolor','none','facecolor', 'cortex'); % plot head model
    alpha 0.9; % opacity of headmodel
    ft_plot_mesh(grid.pos(grid.inside,:)); % plot source model
    grad_mm = ft_convert_units(data.grad, 'mm');
    ft_plot_sens(grad_mm,'style','ob'); % plot MEG channels (sensors)
    hold off;
    view(45,10);

     hf = gcf;
     hf.Position(1:2) = [10 10];
     hf.Position(3:4) = (800 / hf.Position(4)) .* hf.Position(3:4);
     print(hf, [this_output '/step3_source_head_sens_align'], '-dpng', '-r600');
     clear grad_mm
     close all

%%% COMPUTE LEADFIELD -----------------------------------------------------
    % the leadfield is used to provide information on the contribution of a
    % dipole source at a given location in a sensor's region
    cfg              = []; % set up config to prepare the leadfield
    cfg.headmodel    = hdm;
    cfg.sourcemodel  = grid;
    cfg.reducerank   = 2;
    cfg.grad         = data.grad;
    cfg.normalize    = config.step3.normLeadField;
    leadfield        = ft_prepare_leadfield(cfg, data); % create leadfield

%% ACTUAL BEAMFORMING
%%% VECTOR - Time Domain Source Reconstruction ----------------------------

    %%% compute common spatial filter (returns: COVARIANCE MATRIX)
    % the covariance matrix tells us how related the sensors are
    cfg                    = []; % set up config to compute covariance matrix
    cfg.covariance         = 'yes';
    cfg.keeptrials         = 'yes';
    tlock                  = ft_timelockanalysis(cfg, data); % compute covariance matrix

    %%% calculate sensor weights (actual beamforming)
    cfg                 = []; % set up config for sensor weight calculation
    cfg.grad            = data.grad;            % sensor position (gradiometer)
    cfg.headmodel       = hdm;
    cfg.sourcemodel     = grid;          % source model
    cfg.sourcemodel.leadfield  = leadfield.leadfield;
    cfg.lcmv.keepfilter = 'yes';
    source_t_avg        = ft_sourceanalysis(cfg, tlock);
    
    %%% project all trials thru spatial filter
    cfg                  = []; % set up config for beamforming
    cfg.sourcemodel      = grid; % source model
    cfg.sourcemodel.filter = source_t_avg.avg.filter;
    cfg.sourcemodel.leadfield   = leadfield.leadfield;
    cfg.grad             = data.grad; % sensor position (gradiometer)
    cfg.headmodel        = hdm;
    cfg.method           = 'lcmv';
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
        save([ssSubjPath(ss) '/atlas_beamforming_results.mat'],'catmatrix', 'var_explained', 'srate','coords','-mat','-v7.3')
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
