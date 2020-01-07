function [source] = fcp_4_beamforming(paths)
% function [source,source_conn,parc_conn] = fcp_4_beamforming(paths)

%icacleaned data ft_meg_data_cfg.mat
%mri subj/MRI/subj_V2.mri


%3-dimensional source-reconstructed data (not freesurfer - on surface)
% NOTE: the path to the template file is user-specific
%load(fullfile(ftpath, 'template/sourcemodel/standard_sourcemodel3d10mm');
%template_grid = sourcemodel;
%clear sourcemodel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allpath = strsplit(path,':');
% ft_inds = cellfun(@(q) ~isempty(q), strfind(allpath,'fieldtrip'));
% justft = allpath(ft_ind);
% [~,ft_dirind] = min(cellfun(@(x) length(x),justft));
% ftpath = justft(ft_dirind);

config = load_config(paths, paths.name);
config = config.config;
step = 'fcp4';
subj_ds = load_participants(paths,step);
pids = readtable(paths.all_subj_pids);
[subj_match, failure] = ds_pid_match(paths,step);
ssSubjPath = @(x) paths.(subj_match.pid{x});


% NOTE: the path to the template file is user-specific
ftpath = which('ft_defaults.m');
ftpath = ftpath(1:35);
template = ft_read_mri(fullfile(ftpath, '/external/spm8/templates/T1.nii'));
template.coordsys = config.beamforming.template.coordsys;

% segment the template brain and construct a volume conduction model (i.e. head model):
% this is needed to describe the boundary that define which dipole locations are 'inside' the brain.
cfg          = [];
template_seg = ft_volumesegment(cfg, template); % unit: mm

cfg          = [];
cfg.method   = config.beamforming.headmodel.method; %p.sourcemodel.headmodel_method 
template_headmodel = ft_prepare_headmodel(cfg, template_seg); % unit: mm
template_headmodel = ft_convert_units(template_headmodel, config.beamforming.headmodel.units); % Convert the vol to cm, because the CTF convenction is to express everything in cm.

% construct the dipole grid in the template brain coordinates
% the negative inwardshift means an outward shift of the brain surface for inside/outside detection
cfg = [];
cfg.resolution      = config.beamforming.template.grid.resolution;
cfg.tight           = config.beamforming.template.grid.tight;
cfg.inwardshift     = config.beamforming.template.grid.inwardshift;
cfg.headmodel       = template_headmodel;
template_grid       = ft_prepare_sourcemodel(cfg);

% make a figure with the template head model and dipole grid
% figure
% hold on
% ft_plot_vol(template_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
% ft_plot_mesh(template_grid.pos(template_grid.inside,:));
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loas atlas and create a binary mask
atlas=ft_read_atlas([ftpath '/' config.beamforming.atlas.filepath]);
atlas = ft_convert_units(atlas, config.beamforming.headmodel.units);% assure that atlas and template_grid are expressed in the %same units

cfg = [];
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel;
cfg.inputcoord = config.beamforming.atlas.inputcoord;
mask = ft_volumelookup(cfg,template_grid);

% create temporary mask according to the atlas entries
tmp                  = repmat(template_grid.inside,1,1);
tmp(tmp==1)          = 0;
tmp(mask)            = 1;

% define inside locations according to the atlas based mask
template_grid.inside = tmp;

% plot the atlas based grid
% figure;ft_plot_mesh(template_grid.pos(template_grid.inside,:));



rangeOFsubj = 1:length(subj_match.ds);

for ss = rangeOFsubj
    fprintf('\nWorking on subject %s! \n', subj_match.pid{ss});

    % read the single subject anatomical MRI, this should be aligned to MEG head coordinates
    % if the MRI is not aligned, you should use ft_volumerealign
    mri = ft_read_mri([paths.rawmri '/' subj_match.pid{ss} '_V2.mri']);
    mri = ft_convert_units(mri,'cm');
    if any(mri.hdr.fiducial.mri.nas) == 0 || any(mri.hdr.fiducial.mri.lpa) == 0  || any(mri.hdr.fiducial.mri.rpa) == 0
        error('No fiducials found for subject %s!',inddata.ID);
    end

    % segment the anatomical MRI
    cfg        = [];
    cfg.output = 'brain';
    seg        = ft_volumesegment(cfg, mri);

    % construct the volume conductor model (i.e. head model) for each subject
    % this is optional, and for the purpose of this tutorial only required for
    % plotting, later on
    cfg        = [];
    cfg.method = config.beamforming.headmodel.method;
    headmodel_units = 'cm';
    hdm  = ft_prepare_headmodel(cfg, seg);
    hdm = ft_convert_units(hdm,headmodel_units);

    % check segmented mri volumes for proper alignment (brain) and save as image
    seg.transform  = mri.transform;
    seg.anatomy    = mri.anatomy;

    cfg = [];
    cfg.method          = config.beamforming.checkMRIvolumes.method;
    cfg.slicesdim       = config.beamforming.checkMRIvolumes.slicesdim;
    cfg.nslices         = config.beamforming.checkMRIvolumes.nslices;
    % cfg.anaparameter    = 0.5;
    % cfg.title           = ['Segmentation: ', inddata.ID];
%     ft_sourceplot(cfg, seg);
% 
%     hf = gcf;
%     hf.Position(1:2) = [10 10]; 
%     hf.Position(3:4) = (800 / hf.Position(4)) .* hf.Position(3:4);
%     print(hf, inddata.segmentationfigure, '-dpng', '-r600');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot of tissue probability maps 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create the subject specific grid, using the template grid that has just been created
    cfg                = [];
    cfg.grid.warpmni   = config.beamforming.subj.grid.warpmni;
    cfg.grid.template  = template_grid;
    cfg.grid.nonlinear = config.beamforming.subj.grid.nonlinear;
    cfg.mri            = mri;
    cfg.grid.unit      =config.beamforming.subj.grid.unit; %mm
    grid               = ft_prepare_sourcemodel(cfg);

    %coords = sourcemodel.pos;

    % make a figure of the single subject headmodel, and grid positions
%     figure; hold on;
%     ft_plot_vol(hdm, 'edgecolor', 'none', 'facealpha', 0.4);
%     ft_plot_mesh(grid.pos(grid.inside,:));

    %load data = ft_meg_data_cfg.mat
    load([ssSubjPath(ss) '/ft_meg_fullyProcessed.mat'],'-mat','data');

%     % check alignment of source model and head model and save as image
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

    %saveas(hf, char(inddata.sourcefigure), 'fig');
    % print(hf, inddata.sourcefigure, '-dpng', '-r600');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% compute the leadfield
    cfg             = [];
    % cfg.grid        = grid; %sourcemodel
    cfg.headmodel   = hdm;
    cfg.channel     = 'MEG';
    cfg.reducerank   = 2;
    cfg.grid.pos     = grid.pos;
    cfg.grid.inside  = grid.inside;
    cfg.normalize    = config.beamforming.leadfield.normalize;
    cfg.grad         = data.grad;
    leadfield        = ft_prepare_leadfield(cfg, data);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% VECTOR - Time Domain Source Reconstruction %%
    % compute common spatial filter (trial covariance matrix)
    cfg = [];
    cfg.covariance         = config.beamforming.timeDomain.covariance;
    cfg.channel            = 'MEG';
    cfg.covariancewindow   = config.beamforming.timeDomain.covariancewindow;
    cfg.vartrllength       = config.beamforming.timeDomain.vartrllength;
    cfg.keeptrials         = config.beamforming.options.keeptrials;
    tlock = ft_timelockanalysis(cfg, data);  % calculates covariance matrix
    % calculate sensor weights (actual beamforming)
    cfg = [];
    cfg.grad            = data.grad;            % sensor position (gradiometer)
    cfg.headmodel       = hdm;
    cfg.grid            = grid;          % source model
    cfg.grid.leadfield  = leadfield.leadfield;
    cfg.keepfilter      = config.beamforming.options.keepfilter;
    source_t_avg = ft_sourceanalysis(cfg, tlock);
    
    % project all trials through common spatial filter
    cfg = [];
    cfg.grid             = grid; % source model
    cfg.grid.filter      = source_t_avg.avg.filter;
    cfg.grid.leadfield   = leadfield.leadfield;
    cfg.grad             = data.grad; % sensor position (gradiometer)
    %cfg.vol              = hdm; % head model
    cfg.headmodel = hdm;
    cfg.method           = config.beamforming.method;
    cfg.keeptrials       = config.beamforming.options.keeptrials;
    cfg.keepfilter       = config.beamforming.options.keepfilter;
    cfg.rawtrial         = config.beamforming.options.rawtrial;
    source_t_trials = ft_sourceanalysis(cfg, tlock);

    % project the virtual sources to their strongest (dominant) orientation
    % (taking the largest eigenvector of the sources timeseries)
    cfg = [];
    cfg.projectmom = config.beamforming.timeDomain.projectmom;
    cfg.keeptrials = 'yes';
    projection = ft_sourcedescriptives(cfg, source_t_trials);
    
    %%% Interpolate AAL atlas onto Virtual Sources - know which sources are part of which AAL %%%
    % Setup for AAL interpolation
    sourcemodel.pos=template_grid.pos; % allows for interpolation and connectivity analysis

    fullPath = which('ft_preprocessing.m');
    [pathstr,~,~] = fileparts(fullPath);
    atlas=ft_read_atlas([pathstr,'/template/atlas/aal/ROI_MNI_V4.nii']);
    atlas=ft_convert_units(atlas,'cm');

    cfg=[];
    cfg.interpmethod='nearest';
    cfg.parameter='tissue';
    source_aal=ft_sourceinterpolate(cfg,atlas,sourcemodel);

    
    vector_virtual_sources=[];
    ori_avg=[];
    for t=1:projection.df
        vector_virtual_sources_trial=[];
        ori_avg_trial=[];
        ori_116=[];
        for i=1:90 % for each AAL region
            source_timeseries=[];
            ori_region=[];
            aal_node=find(source_aal.tissue==i);
            sources_innode=projection.trial(t).mom(aal_node);
            % creates average time series
            ori_innode=projection.trial(t).ori(aal_node);
            for j=1:length(sources_innode)
                source_timeseries=cat(1,source_timeseries,sources_innode{j}); % matrix of all source time series in region
                ori_region=cat(1,ori_region,ori_innode{1,j});
            end
            if length(sources_innode)>1
                vector_virtual_sources_trial=cat(1,vector_virtual_sources_trial,mean(source_timeseries));
                ori_avg_trial=cat(1,ori_avg_trial,mean(ori_region));
            elseif isempty(sources_innode)
                error('NO NODE %d\n',i);
            else
                vector_virtual_sources_trial=cat(1,vector_virtual_sources_trial,source_timeseries);
                ori_avg_trial=cat(1,ori_avg_trial,ori_region);
            end
            ori_region_padded=zeros(ceil(0.03*nnz(source_aal.tissue)),3); % can have an error if resolution of sourcemodel changed
            ori_region_padded(1:size(ori_region,1),:)=ori_region;
            ori_116=cat(3,ori_116,ori_region_padded); % orientation of eigenvectors of each source (divided into AAL regions)
        end
        vector_virtual_sources=cat(3,vector_virtual_sources,vector_virtual_sources_trial);
        ori_avg=cat(3,ori_avg,ori_avg_trial); % orientation of eigenvector in each AAL region (same for all trials)
    end
    catmatrix=permute(vector_virtual_sources,[2 3 1]); % to have exact same dimension order as old preprocessing (catmatrix)
    srate=data.fsample;
    coords = projection.pos;
    save([ssSubjPath(ss) '/AAL_beamforming_results'],'catmatrix','srate','coords','-mat','-v7.3')
        
    clear coords catmatrix srate vector_virtual_sources vector_virtual_sources_trial source_timeseries ...
        sources_innode ori_avg ori_avg_trial ori_116 atlas sourcemodel source_t_trials projection seg mri data grid hdm ...
        source_t_avg tlock leadfield aal_node ori_region ori_innode ori_region_padded
    fprintf('\nDone subject %s! \n', subj_match.pid{ss});
end
