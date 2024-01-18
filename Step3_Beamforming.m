function Step3_Beamforming(config, pid, visitnum)

% to fill this part in...

%% SET UP LOGGING FILE

if exist('visitnum', 'var')
    this_output = [config.meta.project_path '/' config.meta.analysis_name '/' pid '/' sprintf('ses-%.2d', visitnum)]; % indicate subject-specific output folder path
else
    this_output = [config.meta.project_path '/' config.meta.analysis_name '/' pid]; 
end

load([this_output '/out_struct.mat'])

%% TEMPLATE SOURCE MODEL

% load from FieldTrip templates
load(sprintf('%s/template/sourcemodel/standard_sourcemodel3d%dmm', ...
                config.meta.fieldtrip_path, ...
                config.step3.templateRes), 'sourcemodel');
template_grid = ft_convert_units(sourcemodel, 'mm');
clear sourcemodel;

%% PARTICIPANT MODEL

%%% LOAD ANATOMICAL MRI DATA ----------------------------------------------
    T1_path = [config.meta.project_path '/Prep_T1s/' pid '/' sprintf('ses-%.2d', visitnum)];
    load([T1_path '/Prep_T1_aligned.mat']);
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

 %% INTERPOLATE ATLAS ONTO VIRTUAL SOURCES ----------------------------

    % setup for atlas interpolation - get coordinates
    if config.step3.nativeSpace % when atlas MRI file is already native
        sourcemodel.pos = grid.pos;
    else
        sourcemodel.pos = template_grid.pos; % MNI space atlas
    end
    
    % load atlas
    labeltype = 'tissue';
    if contains(config.step3.atlas, 'mmp') % if MMP glasser atlas in MNI
        atlas = ft_read_atlas([config.meta.megneto_path '/external/atlas/mmp.mat']);
    elseif contains(config.step3.atlas, 'aal')
        atlas                           = ft_read_atlas([config.meta.fieldtrip_path '/template/atlas/aal/ROI_MNI_V4.nii']);
        atlas.tissuelabel               = atlas.tissuelabel(1:90); % we only want non-cerebellar regions (isolate desired regions)
        atlas.tissue(atlas.tissue > 90) = 0;
    elseif contains(config.step3.atlas, 'wmp') 
        % MMP native space atlas load
        atlas = ft_read_atlas([config.meta.rawdata_path ...
                                '/derivatives/T1w_fastsurfer_jtseng/' ...
                                pid '/' sprintf('ses-%.2d', visitnum) ...
                                '/mri/hcpmmp1_ordered.mgz']);
        atlas.transformorig     = atlas.transform;
        atlas.transform         = mri.transform;
        atlas.coordsys          = 'ctf';
        atlas_labels            = load([config.meta.megneto_path '/external/atlas/mmp_labels.mat']); % mat file of labels as cellstring
        atlas.parcellationlabel = atlas_labels.hcpmmp1_labels;
        clear atlas_labels
        
        % fMRI-based ROI spehres for L/R-FFA/LOC
        fmri_spheres            = ft_read_atlas([config.meta.project_path ...
                                      '/fMRI_FFA-LOC_ROIs/' ...
                                      pid '_' sprintf('ses-%.2d', visitnum) ...
                                      '_fMRI_spheres.nii.gz']);
        fmri_spheres.transformorig = fmri_spheres.transform;
        fmri_spheres.transform      = mri.transform;
        fmri_spheres.coordsys       = 'ctf';
        fmri_spheres.parcellationlabel = {'L_FFA', 'R_FFA', 'L_LOC', 'R_LOC'}';
        labeltype                   = 'parcellation';
    end

    % source interpolate
    cfg              = [];
    cfg.interpmethod = 'nearest';
    cfg.parameter    = labeltype;
    source_atlas     = ft_sourceinterpolate(cfg,atlas,sourcemodel); % interpolate source activity onto voxels of anatomical description of the brain
    
    % in the case of additional custom sphere ROIs derived from fMRI
    if contains(config.step3.atlas, 'wmp')
        % identify dipoles corresponding to fMRI spheres
        source_fmri = ft_sourceinterpolate(cfg, fmri_spheres, sourcemodel);
    end
    
    % identify pos corresponding to actual ROIs
    if ~(isstring(config.step3.ROIs)) % if it's numbers and not "all"
        roi_pos = any(source_atlas.parcellation == config.step3.ROIs, 2);
        source_atlas.parcellationlabel_all = source_atlas.parcellationlabel;
        source_atlas.parcellationlabel = source_atlas.parcellationlabel(config.step3.ROIs);
        num_rois = length(config.step3.ROIs);
    else % else it says "all" and we should choose any ROI
        roi_pos = any(source_atlas.parcellation > 0, 2);
        num_rois = length(atlas.parcellationlabel);
    end
    
    if contains(config.step3.atlas, 'wmp')
        roi_pos = any([roi_pos, ...
                   any(source_fmri.parcellation > 0, 2)],2);
        
       % handle overlapping dipoles between Glasser ROIs and fMRI
       overlap = find(sum([roi_pos, source_fmri.parcellation] > 0, 2) == 2);
       source_atlas.parcellation(overlap) = 0;
    end

    
    
    % interim notes:
    %   grid = subject specific coordinates
    %   template_grid = template loaded in from fieldtrip
    %   sourcemodel = either grid/template_grid depending on atlas
    %   roi_pos = x,y,z coordinates of dipoles that actually belong to an
    %             ROI: either all ROIs, or a subset of them by index
    
%% COMPUTE LEADFIELD -----------------------------------------------------
    % the leadfield is used to provide information on the contribution of a
    % dipole source at a given location in a sensor's region
    cfg              = []; % set up config to prepare the leadfield
    cfg.headmodel    = hdm;
    cfg.sourcemodel.pos  = grid.pos;
    cfg.sourcemodel.inside = roi_pos;
    cfg.reducerank   = 2;
    cfg.grad         = data.grad;
    cfg.normalize    = config.step3.normLeadfield;
    leadfield        = ft_prepare_leadfield(cfg, data); % create leadfield

%% ACTUAL BEAMFORMING
%%% VECTOR - Time Domain Source Reconstruction ----------------------------
    
    % need to make sure that only MEG-proper channels are selected
    selchan = ft_channelselection({'all', '-MMSTC*'}, data.label);

    %%% compute common spatial filter (returns: COVARIANCE MATRIX)
    % the covariance matrix tells us how related the sensors are
    cfg                    = []; % set up config to compute covariance matrix
    cfg.covariance         = 'yes';
    cfg.keeptrials         = 'yes';
    cfg.channel            = selchan;
    tlock                  = ft_timelockanalysis(cfg, data); % compute covariance matrix

    %%% calculate sensor weights (actual beamforming)
    cfg                 = []; % set up config for sensor weight calculation
    cfg.grad            = data.grad;            % sensor position (gradiometer)
    cfg.headmodel       = hdm;
    cfg.sourcemodel.pos     = sourcemodel.pos;          % source model
    cfg.sourcemodel.inside = leadfield.inside;
    cfg.sourcemodel.leadfield  = leadfield.leadfield;
    cfg.lcmv.keepfilter = 'yes';
    source_t_avg        = ft_sourceanalysis(cfg, tlock);
    
    %%% project all trials thru spatial filter
    cfg                  = []; % set up config for beamforming
    cfg.sourcemodel.pos      = sourcemodel.pos; % source model
    cfg.sourcemodel.inside   = leadfield.inside;
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
    
%% collapse dipoles within ROIs and save data_roi to file

    % set it up in a fieldtrip-looking structure
    data_roi                        = [];
    data_roi.time                   = data.time;
    data_roi.fsample                = data.fsample;
    data_roi.trialinfo              = data.trialinfo;
    data_roi.sourceinterp.pos       = source_atlas.pos;
    data_roi.sourceinterp.(labeltype) = source_atlas.(labeltype);
    data_roi.label                  = source_atlas.parcellationlabel;
    
    for t = 1:projection.df
        %%% AND FOR EACH NODE ---------------------------------------------
        for i = 1:num_rois
            % identify source coords that fall within ROI
            node                     = find(source_atlas.(labeltype)==config.step3.ROIs(i)); 
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
                warning('No sources in ROI %s.\n',source_atlas.([labeltype 'label']){i});
            end
        end
    end
    
    if contains(config.step3.atlas, 'wmp')
        data_roi.label = [data_roi.label; source_fmri.parcellationlabel];
        data_roi.sourceinterp.parcellation(:,2) = source_fmri.parcellation;
        for t = 1:projection.df
        %%% AND FOR EACH NODE ---------------------------------------------
            for i = 1:4
                % identify source coords that fall within ROI
                node                     = find(source_fmri.(labeltype)==i); 
                source_timeseries        = cell2mat(projection.trial(t).mom(node)); % get the timeseries; num_nodes x time
                % ori_region               = cell2mat(projection.trial(t).ori(node)); % orientations; num_nodes x time

                % IF NODE EXISTS
                if size(source_timeseries, 1) >= 1 
                    if config.step3.combineDipoles == "mean"
                        data_roi.trial{t}(num_rois+i,:) = nanmean(source_timeseries,1); % take avg across source points
                    elseif config.step3.combineDipoles == "pca"
                        [~, score, ~, ~, explained] = pca(transpose(source_timeseries)); % perform pca
                        data_roi.trial{t}(num_rois+i,:) = transpose(score(:, 1)); % store first principal component across timeseries
                        data_roi.var{t}(i+1) = explained(1);
                    end
                    % ori_avg(:,t,i) = nanmean(ori_region,1);
                % IF NO SOURCE POINTS W/IN NODE
                else
                    warning('No sources in ROI %s.\n',source_fmri.([labeltype 'label']){i});
                end
            end
        end
    end
    
    save([this_output '/step3_data_roi.mat'], 'data_roi') 
                
end
