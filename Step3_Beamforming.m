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
    % mri     = ft_convert_units(mri,'cm');
    
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
    cfg                 = []; % set up config for participant source model preparation 
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
    cfg.normalize    = config.step3.normLeadfield;
    leadfield        = ft_prepare_leadfield(cfg, data); % create leadfield

%% ACTUAL BEAMFORMING
%%% VECTOR - Time Domain Source Reconstruction ----------------------------
    
    % need to make sure that only MEG-proper channels are selected

    %%% compute common spatial filter (returns: COVARIANCE MATRIX)
    % the covariance matrix tells us how related the sensors are
    cfg                    = []; % set up config to compute covariance matrix
    cfg.covariance         = 'yes';
    cfg.keeptrials         = 'yes';
    cfg.channel           = data.label(1:183);
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
    cfg.projectmom       = 'yes';
    cfg.keeptrials       = 'yes';
    projection           = ft_sourcedescriptives(cfg, source_t_trials); % project to dominant orientation
    
%%% INTERPOLATE ATLAS ONTO VIRTUAL SOURCES ----------------------------

    % setup for atlas interpolation - get coordinates
    sourcemodel.pos = template_grid.pos; 
    
    % load atlas
    if contains(config.step3.atlas, 'mmp') % if MMP glasser atlas
        atlas = ft_read_atlas([config.meta.megneto_path '/external/atlas/mmp.mat']);
    elseif contains(config.step3.atlas, 'aal')
        atlas = ft_read_atlas([config.meta.fieldtrip_path '/template/atlas/aal/ROI_MNI_V4.nii']);
        atlas.tissuelabel   = atlas.tissuelabel(1:90); % we only want non-cerebellar regions (isolate desired regions)
        atlas.tissue(atlas.tissue > 90) = 0;
    end

    % source interpolate
    cfg              = [];
    cfg.interpmethod = 'nearest';
    cfg.parameter    = 'tissue';
    source_atlas       = ft_sourceinterpolate(cfg,atlas,sourcemodel); % interpolate source activity onto voxels of anatomical description of the brain

    % set it up in a fieldtrip-looking structure
    data_roi                        = [];
    data_roi.time                   = data.time;
    data_roi.fsample                = data.fsample;
    data_roi.label                  = atlas.tissuelabel;
    data_roi.sourceinterp.pos       = source_atlas.pos;
    data_roi.sourceinterp.tissue    = source_atlas.tissue;
        
    for t = 1:projection.df
        %%% AND FOR EACH NODE ---------------------------------------------
        for i = 1:max(size(atlas.tissuelabel))
            % identify source coords that fall within ROI
            node                     = find(source_atlas.tissue==i); 
            source_timeseries        = cell2mat(projection.trial(t).mom(node)); % get the timeseries; num_nodes x time
            % ori_region               = cell2mat(projection.trial(t).ori(node)); % orientations; num_nodes x time
            
            % IF NODE EXISTS
            if size(source_timeseries, 1) >= 1 
                if config.step3.combineDipoles == "mean"
                    data_roi.trial{t}(i,:) = nanmean(source_timeseries,1); % take avg across source points
                elseif config.step3.combineDipoles == "pca"
                    [~, score, ~, ~, explained] = pca(transpose(source_timeseries)); % perform pca
                    data_roi.trial{t}(i,:) = transpose(score(:, 1)); % store first principal component across timeseries
                    data_roi.var{t}(i) = explained(1);
                end
                % ori_avg(:,t,i) = nanmean(ori_region,1);
            % IF NO SOURCE POINTS W/IN NODE
            else
                warning('No sources in ROI %s.\n',atlas.tissuelabel{i});
            end
        end
    end
    
    save([this_output '/step3_data_roi.mat'], 'data_roi') 
                
end
