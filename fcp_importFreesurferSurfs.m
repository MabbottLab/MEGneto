function fcp_importFreesurferSurfs(pathsS, paths)

% fcp_importFreesuferSurfs registers CTF mri files with already processed
% freesurfer surf structures and hcp-workbench processing - ready for beamforming 
% 
% NOTES:
%   - Ensure that subj_fcp1.csv is populated with the subject IDs of
%   included participants. 
%   - Desired parameters should be defined in the JSON config file.
%   User should double-check that the JSON config file is populated
%   appropriately, especially if a template JSON was copied over. 
%
% INPUTS:
%   paths               =   struct defining paths to data, participant
%                           folders, analysis folders, config files, etc. 
%
% OUTPUTS:
%
% See also: hcp-workbench,

% Last updated by: Sonya Bells, 2020-02-14
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.

%% SETUP: LOAD CONFIG, PARTICIPANTS, CHECK FOR FULL DATASET, OUTPUTS
% load config JSON with analysis parameters
config      = load_config(paths, paths.name);
config      = config.config;
%%%%%%%%%%%%%%%%step        = 'fcp4';

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

rangeOFsubj = 1:length(subj_match.ds);

for ss = rangeOFsubj

	mripath = fullfile(pathsS.basepath,subj_match.pid{ss});
	datapath = fullfile(mripath,'freesurfer/workbench');
	freesurferpath = fullfile(mripath,'freesurfer/mri');

	%%% FOR EACH PARTICIPANT --------------------------------------------------
    fprintf('\nWorking on subject %s! \n', subj_match.pid{ss});

	%%% LOAD ctf mri image - fiducials marked ---------------------------------

	%mri_orig = ft_read_mri([sub,'/',sub'_T1_V2.mri']);
	mri = ft_read_mri([paths.rawmri '/' subj_match.pid{ss} '_T1_V2.mri'])
	% mri = ft_determine_coordsys(mri_orig, 'interactive', 'yes');
	% mri_orig = ft_determine_coordsys(mri_orig, 'interactive', 'no');
    
    % check for fiducials
    if any(mri.hdr.fiducial.mri.nas) == 0 || any(mri.hdr.fiducial.mri.lpa) == 0  || any(mri.hdr.fiducial.mri.rpa) == 0
        error('No fiducials found for subject %s!', subj_match.pid{ss});
    end

	%%% LOAD FreeSurfer T1 ----------------------------------------------------
	t1 = ft_read_mri([freesurferpath,'/norm.mgz']);
	t1 = ft_determine_coordsys(t1,'interactive', 'no');
	%rasa
	

	%%% Register FreeSurfer T1 -> CTR MRI (sensor space) -----------------------
	cfg                     = [];
	cfg.method        = 'spm';
	cfg.spmversion  = 'spm12';
	cfg.coordsys      = 'ctf';
	cfg.viewresult    = 'yes';
	%ft_volumerealign(cfg, mri, target)
	mri_t12ctf = ft_volumerealign(cfg, t1, mri);

	%save transform
	transform_t12ctf = mri_t12ctf.transform; 
	save(fullfile(mripath,sprintf('%s_transform_t12ctf_n',subj_match.pid{ss})), 'transform_t12ctf');


	%save registered volume
	%%
	cfg          = [];
	cfg.filename = fullfile(mripath,sprintf('%s_t12ctf_n', subj_match.pid{ss}));
	cfg.filetype = 'nifti';
	cfg.parameter = 'anatomy';
	ft_volumewrite(cfg, mri_t12ctf);

	transform_surf2mri = mri_t12ctf.transform/mri_t12ctf.transformorig;

	%%% Load in FreeSurfer surfaces -----------------------
	filename = fullfile(datapath,[subj_match.pid{ss},'.L.midthickness.8k_fs_LR.surf.gii']);
	sourcemodel_8korig = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});


	sourcemodel_8k = ft_transform_geometry(transform_surf2mri, sourcemodel_8korig);
	sourcemodel_8k.inside = sourcemodel_8k.atlasroi>0;
	sourcemodel_8k = rmfield(sourcemodel_8k, 'atlasroi');
	
	%%% Create Headmodel -------------------------------

	%load template_grid (can create your own template grid)
	% NOTE: the path to the template file is user-specific
	a = load(fullfile(pathsS.ftpath, 'template/sourcemodel/standard_sourcemodel3d10mm'));
	template_grid = a.sourcemodel;
	% %clear sourcemodel;

	%%
	% segment the anatomical CTF MRI
	cfg        = [];
	cfg.output = 'brain';
	seg        = ft_volumesegment(cfg, mri);


	%%
	% construct the volume conductor model (i.e. head model) for each subject	
	% this is optional, and for the purpose of this tutorial only required for
	% plotting, later on
	cfg        = [];
	cfg.method = 'singleshell';
	headmodel  = ft_prepare_headmodel(cfg, seg);
		
	%%
	cfg = [];
	cfg.intersectmesh = headmodel.bnd;
	ft_sourceplot(cfg, mri);

	%%
	% create the subject specific grid, using the template grid that has just been created
	cfg           = [];
	cfg.warpmni   = 'yes';
	cfg.template  = template_grid;
	cfg.nonlinear = 'yes';
	cfg.mri       = mri;
	cfg.unit      ='mm';
	grid          = ft_prepare_sourcemodel(cfg);

	%%% Save Headmodel and surface 8k -------------------------------
	save(fullfile(mripath,sprintf('%s_sourcemodel_8k',subj_match.pid{ss})), 'sourcemodel_8k');
	save(fullfile(mripath,sprintf('%s_hdm',subj_match.pid{ss})), 'headmodel');

	%%% LOAD MEG DATA ---------------------------------------------------------
    load([ssSubjPath(ss) '/ft_meg_fullyProcessed.mat'],'-mat','data');

	figure;
	% make the headmodel surface transparent
	ft_plot_headmodel(headmodel, 'edgecolor', 'none'); alpha 0.9
	ft_plot_mesh(ft_convert_units(sourcemodel_8k, 'mm'),'vertexcolor',sourcemodel_8k.sulc);
	ft_plot_sens(data.grad);
	view([0 -90 0])

	%%% Prepare leadfield -------------------------------
	

	 %% compute the leadfield
	cfg 			= [];
	cfg.grid 		= sourcemodel;
	cfg.headmodel 	= headmodel;
	cfg.channel 	= {'MEG'};
	lf 				= ft_prepare_leadfield(cfg, data);



	%%% Source reconstructrion -------------------------------
	
	%%% ? -------------------------------
	
	%%%?  -------------------------------
	
	% In addistion to forward model, the beamforming needs a sensor-level covariance matrix
	% or a cross-spectral density matrix

	% Compute sensor level Fourier spectra, to be used for cross-spectral density computation
	cfg 			= [];
	cfg.method		= 'mtmfft'; % analysis an entire spectrum for entire data length (mulitaper freq transformation)
	cfg.output 		= 'fourier'; % return complex Fourier-spcetra ; powandcsd return the power and cross-spectra
	%cfg.taper = 'dpss'; %default 
	cfg.keeptrials 	= 'yes';
	cfg.tapsmofrq 	=1;
	cfg.foi 		=10;
	freq 			= ft_freqanalysis(cfg,data)

% cfg.method = 'mtmconvol'; %implements mulitaper time-frequency transformation based on multiplication in freq domain
% cfg.taper  = 'hanning';
% cfg.foi 	 = 1:1:50; % analysis 1 to 50Hz in step sof 1Hz (frequencies of interest)
% cfg.t_ftimwin = ones(length(cfg.foi),1).*0.03; %length of time window in seconds = 0.03s
% cfg.toi  	= 0:0.05:6.99; %the times of which the analysis windows should be centered (in seconds) from 0 to 6.99s in steps of 0.05
% 5./cfg.foi
% cfg.output = 'powandcsd'; %return power and cross-spectra
% cfg.tapsmofrq = 2; % number, the amount of spectral smoothing through muli-tapering (e.g. 2Hz) (0.4*cfg.foi)
% cfg.keeptrials = 'yes';
% freq_wpli = ft_freqanalysis(cfg, data);


	% Source reconstructrion	
	% Using method 'pcc' Implements DICS (underlying algorithm for computing the spatial filters is according to DICS)
	% , but provides more flexibility with respenct to data handing
	cfg 				= [];
	cfg.frequency 		= freq.freq;
	cfg.method 			= 'pcc'; %DICS and PCC for freq or time-freq data
	%cfg.method 		= 'eloreta' %used both for time, freq and time-freq
	cfg.grid 			= lf;
	cfg.headmodel 		= headmodel;
	cfg.keeptrials 		= 'yes';
	cfg.pcc.lambda 		= '10%';
	cfg.pcc.projectnoise = 'yes';
	cfg.pcc.fixedori 	='yes';
	source = ft_sourceanalysis(cfg, freq)
	source = ft_sourcedescriptives([], source) % to get NAI

	%Visualize NAI (reconstructed activity (NAI) of resting state alpha power)
	cfg 				= [];
	cfg.method 			= 'surface';
	cfg.funparameter 	= 'nai';
	cfg.maskparameter 	= cfg.funparameter;
	cfg.funcolorlim 	= [0.0 8];
	cfg.opacitylim 		= [3 8];
	cfg.opacitymap 		= 'rampup';
	cfg.funcolormap 	='jet';
	cfg.colorbar 		= 'no';
	ft_sourceplot(cfg, source);
	view([-90 30]);
	light;
	
	%%%?  -------------------------------
	%%%?  -------------------------------


end