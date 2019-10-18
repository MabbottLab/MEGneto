function [source,source_conn,parc_conn] = fcp_4_beamforming(paths)

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
ftpath = ftpath(1:30);
template = ft_read_mri(fullfile(ftpath, '/external/spm8/templates/T1.nii'));
template.coordsys = config.beamforming.template.coordsys;
% template.coordsys = config.beamforming.template.coordsys; % so that FieldTrip knows how to interpret the coordinate system

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
figure
hold on
ft_plot_vol(template_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
ft_plot_mesh(template_grid.pos(template_grid.inside,:));

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
figure;ft_plot_mesh(template_grid.pos(template_grid.inside,:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For visualization only

rangeOFsubj = 1:length(subj_match.ds);

for ss = rangeOFsubj
    % read the single subject anatomical MRI, this should be aligned to MEG head coordinates
    % if the MRI is not aligned, you should use ft_volumerealign
    mri = ft_read_mri([paths.rawmri '\' subj_match.pid{ss} '_V2.mri']);
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
    hdm  = ft_prepare_headmodel(cfg, seg);
    hdm = ft_convert_units(hdm,headmodel_units);


    % check segmented mri volumes for proper alignment (brain) and save as image
    seg.transform  = mri.transform;
    seg.anatomy    = mri.anatomy;

    cfg = [];
    cfg.method          = config.beamforming.checkMRIvolumes.method;
    cfg.slicesdim       = config.beamforming.checkMRIvolumes.slicesdim;
    cfg.nslices         = config.beamforming.checkMRIvolumes.nslices;
    cfg.maskparameter    = 0.5;
    % cfg.title           = ['Segmentation: ', inddata.ID];
    ft_sourceplot(cfg, seg);

    hf = gcf;
    hf.Position(1:2) = [10 10]; 
    hf.Position(3:4) = (800 / hf.Position(4)) .* hf.Position(3:4);
    print(hf, inddata.segmentationfigure, '-dpng', '-r600');
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
    figure; hold on;
    ft_plot_vol(hdm, 'edgecolor', 'none', 'facealpha', 0.4);
    ft_plot_mesh(grid.pos(grid.inside,:));

    %load data = ft_meg_data_cfg.mat

    % check alignment of source model and head model and save as image
    figure('Name',['Subject: ', inddata.ID]);
    hold on;
    ft_plot_vol(hdm,'edgecolor','none','facecolor', 'cortex'); % plot head model
    alpha 0.4; % opacity of headmodel
    ft_plot_mesh(grid.pos(grid.inside,:)); % plot source model
    ft_plot_sens(data.grad,'style','ob'); % plot MEG channels (sensors)
    hold off;
    view(45,10);

    hf = gcf;
    hf.Position(1:2) = [10 10];
    hf.Position(3:4) = (800 / hf.Position(4)) .* hf.Position(3:4);

    %saveas(hf, char(inddata.sourcefigure), 'fig');
    print(hf, inddata.sourcefigure, '-dpng', '-r600');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% compute the leadfield
cfg             = [];
cfg.grid        = grid; %sourcemodel
cfg.headmodel   = hdm;
cfg.channel     = {'MEG'};
cfg.reducerank   = config.beamforming.leadfield.reducerank;
cfg.normalize    = config.beamforming.leadfield.normalize;
cfg.grad         = data.grad;
lf              = ft_prepare_leadfield(cfg, data);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(inddata.process,'Tlock') 
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
    cfg.headmodel = hdm;
    cfg.grid            = grid;          % source model
    cfg.grid.leadfield  = lf.leadfield;
    cfg.keepfilter      = config.beamforming.options.keepfilter;
    source_t_avg = ft_sourceanalysis(cfg, tlock);
    
    % project all trials through common spatial filter
    cfg = [];
    cfg.grid             = grid; % source model
    cfg.grid.filter      = source_t_avg.avg.filter;
    cfg.grid.leadfield   = lf.leadfield;
    cfg.grad             = data.grad; % sensor position (gradiometer)
    %cfg.vol              = hdm; % head model
    cfg.headmodel = vol;
    cfg.method           = p.beamformer.method;
    cfg.keeptrials       = config.beamforming.options.keeptrials;
    cfg.rawtrial         = config.beamforming.timeDomain.rawtrial;
    cfg.keepfilter       = config.beamforming.options.keepfilter;
    source_t_trials = ft_sourceanalysis(cfg, tlock);

    % project the virtual sources to their strongest (dominant) orientation
    % (taking the largest eigenvector of the sources timeseries)
    cfg = [];
    cfg.projectmom = config.beamforming.timeDomain.projectmom;
    cfg.keeptrials = 'yes';
    projection = ft_sourcedescriptives(cfg, source_t_trials);

elseif strcmp(inddata.process,'csd')

	%% compute sensor level Fourier spectra, to be used for cross-spectral density computation.
	cfg            = [];
	cfg.method     = config.beamforming.fourier.method;
	cfg.output     = config.beamforming.fourier.output;
	cfg.keeptrials = 'yes';
    cfg.channel     = {'MEG'};
	cfg.tapsmofrq  = config.beamforming.fourier.tapsmofrq;
	cfg.foi      =config.beamforming.fourier.foi; %inddata.freq;
	freq           = ft_freqanalysis(cfg, data);

	%pcc’ as method : implements DICS (the underlying algorithm for computing the spatial filters is according to DICS)
	%the advantage is that the ‘pcc’-implementation directly outputs, for each dipole location in the sourcemodel, 
	%the fourier coefficients (i.e. phase and amplitude estimates) for each of the trials. This can subsequently be used 
	%in a straightforward way for connectivity analysis

	%% do the source reconstruction
	cfg                   = [];
	cfg.frequency         = freq.freq; %for single Hz
	cfg.method            = config.beamforming.sourceReconstruction.method; %partial cannonical correlation/coherence
    cfg.grad              = freq.grad;
	cfg.grid              = lf;
	cfg.headmodel         = hdm;
	cfg.keeptrials        = config.beamforming.options.keeptrials;
	cfg.pcc.lambda        = config.beamforming.sourceReconstruction.pcc.lambda;
	cfg.pcc.projectnoise  = config.beamforming.sourceReconstruction.pcc.projectnoise;
	cfg.pcc.fixedori      = config.beamforming.sourceReconstruction.pcc.fixedori;
	source = ft_sourceanalysis(cfg, freq);
    
    %% reduce the source reconstructed data to the dominant orientation
    cfg = [];
    cfg.projectmom = 'yes';
	source_proj = ft_sourcedescriptives(cfg, source); % to get the neural-activity-index
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%interpolate the functional data onto the anatomical data

% cfg = [];

%mri_resliced = ft_volumereslice(cfg, mri);

%neural activity index (NAI), in order to remove the center of the head bias shown above. 
cfg = [];
cfg.downsample = config.beamforming.interpolate.downsample;
cfg.parameter = config.beamforming.interpolate.parameter;
% cfg.parameter = 'avg.nai';
source_powinterp = ft_sourceinterpolate(cfg, source_proj, mri);

%read atlas
atlas=ft_read_atlas([ftpath '/' config.beamforming.atlas.filepath]);
atlas = ft_convert_units(atlas,'cm');% assure that atlas and template_grid are expressed in the %same units
figure;imagesc(atlas.tissue(:,:,68));

cfg=[];
cfg.interpmethod= config.beamforming.interpolate.method;
cfg.parameter= config.beamforming.atlas.parameter;
source_bna=ft_sourceinterpolate(cfg,atlas,source_proj);




grid.pos=template_grid.pos;
%read atlas
atlas=ft_read_atlas([ftpath '/' config.beamforming.atlas.filepath]);
atlas=ft_convert_units(atlas,'cm');

cfg=[];
cfg.interpmethod= config.beamforming.interpolate.method;
cfg.parameter= config.beamforming.atlas.parameter;
mri_bna=ft_sourceinterpolate(cfg,atlas,grid);


cfg=[];
cfg.interpmethod= config.beamforming.interpolate.method;
cfg.parameter= config.beamforming.atlas.parameter;
source_bna=ft_sourceinterpolate(cfg,atlas,source_proj);
% indx = find(sourcemodel2.tissue==80);


cfg            = [];
cfg.downsample = config.beamforming.interpolate.downsample;
cfg.parameter = config.beamforming.interpolate.parameter;
% cfg.parameter = 'avg';
cfg.interpmethod = config.beamforming.interpolate.method;
sourcePostInt_nocon  = ft_sourceinterpolate(cfg, source_proj , atlas);


%% compute connectivity
cfg         = [];
cfg.method  =inddata.connectivitymethod; %'coh';
cfg.complex = config.beamforming.connectivity.complex;
source_conn2 = ft_connectivityanalysis(cfg, sourcePostInt_nocon);

figure;imagesc(source_conn2.cohspctrm);




cfg            = [];
% cfg.downsample = 2;
cfg.interpmethod = config.beamforming.interpolate.method;
cfg.parameter = config.beamforming.connectivity.parameter;
sourcePostInt  = ft_sourceinterpolate(cfg, source_conn , atlas);







