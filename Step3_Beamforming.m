function Step3_Beamforming(config, pid, visitnum)

% Step3_Beamforming carries out beamforming and source projection on
% cleaned data. Steps:
%       1. Template grid (sourcemodel) from FieldTrip is loaded in with preselected
%           dipole grid resolution, in MNI position space
%       2. Participant model is created using the aligned T1 (output of Prep_T1)
%           to create a brain segmentation for use in modelling
%       3. Image saved to file depicting alignment between head, source, and
%           sensors
%       4. Atlas is loaded in and dipole grid is interpolated to atlas,
%           thereby assigning to each dipole/source an atlas label. Based
%           on user input, only dipoles within relevant ROIs are source
%           reconstructed to save on compute resources
%       5. Beamforming occurs with all sensor weight calculations. 
%       6. Dipoles all belonging to particular ROIs are combined into
%           representative timeseries for that ROI, then saved.
%
% INPUTS:
%   > config: 
%       struct, configured in your "main" script with
%       all analysis parameters and options and paths
%   > pid:
%       string, participant ID used to build the paths
%       to relevant data files and output folders
%   > visitnum:
%       int, visit number if longitudinal, optional arg)
%
% OUTPUTS:
%   > step3_data_roi.mat:
%       fieldtrip style data object with source-reconstructed data by ROI
%   > step3_source_head_sens_align.png:
%       image of the lineup between sensors, dipole sources, brain
%       segmentation
%
% NOTES ON VARIABLES BELOW:
%   grid = subject specific coordinates
%   template_grid = template loaded in from fieldtrip
%   sourcemodel = either grid/template_grid depending on atlas
%   roi_pos = x,y,z coordinates of dipoles that actually belong to an
%             ROI: either all ROIs, or a subset of them by index
%
% See also: FT_READ_MRI, FT_VOLUMESEGMENT, FT_CONVERT_UNITS, FT_PREPARE_HEADMODEL,
% FT_PREPARE_SOURCEMODEL, FT_RESAMPLEDATA, FT_PREPARE_LEADFIELD,
% FT_TIMELOCKANALYSIS, FT_SOURCEANALYSIS, FT_SOURCEDESCRIPTIVES,
% FT_READ_ATLAS, FT_SOURCEINTERPOLATE, FT_VOLUMELOOKUP
%
% Last updated by: Julie Tseng, 2024-09-09
%   This file is part of MEGneto, see https://github.com/MabbottLab/MEGneto
%   for the documentation and details.

%% SETUP

if exist('visitnum', 'var')
    this_output = [config.meta.project_path '/' config.meta.analysis_name '/' ...
        pid '/' sprintf('ses-%.2d', visitnum)]; % indicate subject-specific output folder path
else
    this_output = [config.meta.project_path '/' config.meta.analysis_name '/' pid]; 
end

load([this_output '/out_struct.mat'])

%% TEMPLATE SOURCE MODEL

% load from FieldTrip templates a dipole grid with user-specified
% resolution (defines spacing between dipoles, e.g., 5mm)
load(sprintf('%s/template/sourcemodel/standard_sourcemodel3d%dmm', ...
                config.meta.fieldtrip_path, ...
                config.step3.templateRes), 'sourcemodel');

% convert units and clean up
template_grid = ft_convert_units(sourcemodel, 'mm');
clear sourcemodel;

%% PARTICIPANT MODEL

%%% LOAD ANATOMICAL MRI DATA ----------------------------------------------
    if exist('visitnum', 'var')
        this_T1_path = fullfile(config.step3.PrepT1path, pid, ...
                        sprintf('ses-%.2d', visitnum), "Prep_T1_aligned.mat");
    else
        this_T1_path = fullfile(config.step3.PrepT1path, pid, ...
                "Prep_T1_aligned.mat");
    end
    load(this_T1_path);
    
    % check for fiducials which help to localize head position relative to
    % the sensors
    if ~isfield(mri.cfg.fiducial, 'nas')
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
    elseif contains(config.step3.atlas, 'glasser_native') 
        % MMP native space atlas load
        
        % first, build the path to fastsurfer output
        atlas_path = fullfile(config.step3.glasserPath, ...
                                pid, sprintf('ses-%.2d', visitnum), ...
                                '/mri/hcpmmp1_ordered.mgz');
        
        % check for the file's existence and throw error otherwise
        if ~isfile(atlas_path)
            error("Did not find the hcpmmp1_ordered.mgz native space Glasser parcellation in the glasserPath folder - please fix!")
        end
        
        % if it exists, load it in
        atlas = ft_read_atlas(atlas_path);

        % assuming the native atlas parcellation and your Prep_T1 are from the same space, 
        % pull over the spatial transforms that were already identified
        % during Prep_T1
        atlas.transformorig     = atlas.transform;
        atlas.transform         = mri.transform;
        atlas.coordsys          = 'ctf';
        
        % grab ROI labels
        atlas_labels            = load([config.meta.megneto_path '/external/atlas/mmp_labels.mat']); % mat file of labels as cellstring
        atlas.parcellationlabel = atlas_labels.hcpmmp1_labels;
        clear atlas_labels
    end

    % source interpolate: this takes each dipole grid position and assigns
    % to it an atlas label based on the lineup between the sourcemodel
    % (template grid) and the atlas - critical at this point that the
    % sourcemodel and atlas are both aligned in positions
    cfg              = [];
    cfg.interpmethod = 'nearest';
    cfg.parameter    = labeltype;
    source_atlas     = ft_sourceinterpolate(cfg,sourcemodel,atlas); 
   
    % identify dipoles corresponding to actual ROIs
    if ~(isstring(config.step3.ROIs)) % if it's numbers and not "all"
        % find corresponding labels to specified ROI indices
        roi_labels = atlas.parcellationlabel(config.step3.ROIs);
        
        % ensure index values are the exact same
        roi_idx = find(ismember(source_atlas.parcellationlabel, ...
                                string(atlas.parcellationlabel(config.step3.ROIs))));
        
        % find dipoles that fall into one of the specified ROIs                    
        roi_pos = any(source_atlas.parcellation == roi_idx', 2);
        
        % retain list of all ROIs
        source_atlas.parcellationlabel_all = source_atlas.parcellationlabel;
        
        % grab number of ROIs
        num_rois = length(config.step3.ROIs);
    else % else it says "all" and we should choose any ROI
        roi_pos = any(source_atlas.parcellation > 0, 2);
        num_rois = length(atlas.parcellationlabel);
        config.step3.ROIs = 1:num_rois; % <----- this is the new line to try adding
    end
    
    
%% COMPUTE LEADFIELD -----------------------------------------------------
    % the leadfield is used to provide information on the contribution of a
    % dipole source at a given location in a sensor's region
    cfg                         = []; % set up config to prepare the leadfield
    cfg.headmodel               = hdm;
    cfg.sourcemodel.pos         = grid.pos;
    cfg.sourcemodel.inside      = roi_pos; % defines which sources to reconstruct
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
    cfg                      = []; % set up config for beamforming
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
    
    for t = 1:projection.df % FOR EACH TRIAL-------------------------------
        %%% AND FOR EACH NODE ---------------------------------------------
        for i = 1:num_rois
            % handle filtered ROIs
            this_roi_idx = find(string(source_atlas.(sprintf('%slabel_all', labeltype))) == roi_labels{i});
            data_roi.label(i,1) = string(roi_labels{i});
            
            % identify source coords that fall within ROI
            node                     = find(source_atlas.(labeltype)==this_roi_idx); 
            source_timeseries        = cell2mat(projection.trial(t).mom(node)); % get the timeseries; num_nodes x time
            
            % IF NODE EXISTS
            if size(source_timeseries, 1) >= 1 
                if config.step3.combineDipoles == "mean"
                    data_roi.trial{t}(i,:) = nanmean(source_timeseries,1); % take avg across source points
                elseif config.step3.combineDipoles == "pca"
                    [~, score, ~, ~, explained] = pca(transpose(source_timeseries)); % perform pca
                    data_roi.trial{t}(i,:) = transpose(score(:, 1)); % store first principal component across timeseries
                    data_roi.var{t}(i) = explained(1);
                end
            % IF NO SOURCE POINTS W/IN NODE
            else
                warning('No sources in ROI %s.\n',source_atlas.([labeltype 'label']){this_roi_idx});
                data_roi.trial{t}(i,:) = NaN(1, size(data_roi.time{1}, 2));
            end
        end
    end
    
    save([this_output '/step3_data_roi.mat'], 'data_roi') 
                
end
