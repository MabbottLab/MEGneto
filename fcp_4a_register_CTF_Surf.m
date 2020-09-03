function fcp_4a_register_CTF_Surf(paths, surf_path)

% fcp_4_Beamforming_Freq_Surf does the following:
% 1) registers CTF mri files with already processed
% freesurfer surf structures and hcp-workbench processing
% 2) frequency analysis on cleaned data
% 3) then carries out beamforming and source projection 
%
%
% NOTES:
%   - Ensure that subj_fcp4.csv is populated with the subject IDs of
%   participants you want to include after checking over initial results. 
%   - Assumes that anatomical MRI is already aligned with MEG data. If not
%   aligned, use ft_volumerealign function. 
%
% INPUTS:
%   paths               =   struct defining paths to data, participant
%                           folders, analysis folders, config files, etc. 
%
% OUTPUTS:
%   .mat
%       ?? - catmatrix:    num_samples x num_trials x num_nodes
%       ?? - srate:        sampling rate
%       ?? - coords:       coordinates of source points w/in
%
% See also: 
%
% Created by: Sonya Bells, 2020-02-14
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

% initialize output files
    % images
%         fcp4_output.fig_headsurf  = 'headsurf.png';

%% PARTICIPANT MODELS

rangeOFsubj = 1:length(subj_match.ds);

for ss = rangeOFsubj
    %%% PREPARE FOR COREGISTRATION-----------------------------------------
    fprintf('\nWorking on subject %s! \n', subj_match.pid{ss});

    % set up paths to MRI, freesurfer surface, data locations
    mripath         = fullfile(surf_path,subj_match.pid{ss});
    datapath        = fullfile(mripath,'freesurfer/workbench');
    freesurferpath  = fullfile(mripath,'freesurfer/mri');

    % load ctf image (mm units)
    mri             = ft_read_mri([paths.rawmri '/' subj_match.pid{ss} '_T1_V2.mri']);

    % check for fiducials
    if any(mri.hdr.fiducial.mri.nas) == 0 || any(mri.hdr.fiducial.mri.lpa) == 0  || any(mri.hdr.fiducial.mri.rpa) == 0
        error('No fiducials found for subject %s!', subj_match.pid{ss});
    end

    % mark AC space for better co-registration compatibility
        % step 1: use three images to find tip of anterior commissure
        % step 2: press 'a' on keyboard to mark ac 
        %         press 'p' on keyboard to mark pc 
        %         press 'z' on keyboard to mark superior / up
        %         press 'r' on keyboard to mark right side
        % step 3: press 'f' to check where you marked the point
        %         press 'c' to toggle the crosshair
        % step 4: press 'q' to finalize the markers and quit
    cfg             = [];
    cfg.method      = 'interactive';
    cfg.coordsys    = 'acpc';
    cfg.spmversion  = 'spm12';
    mri_acpc        = ft_volumerealign(cfg,mri);ac 

	% load FreeSurfer MRI, then specify coordinate system
    % need to indicate ras coordinate system (right, ant, superior)
    t1              = ft_read_mri([freesurferpath, '/norm.mgz']);
    t1              = ft_determine_coordsys(t1, 'interactive', 'yes');
 
%     % skull-strip  the anatomical CTF MRI
%     cfg               = [];
%     cfg.output        = 'skullstrip';
%     cfg.viewresult    = 'yes';
%     mri_acpc_brain    = ft_volumesegment(cfg, mri_acpc);
   
    %%% REGISTER FreeSurfer T1 -> CTR MRI ---------------------------------
    cfg                 = [];
    cfg.method          = 'spm';
    % cfg.spmversion = 'spm8';
    cfg.spmversion      = 'spm12';
    cfg.viewresult      = 'yes';
    mri_t12acpc         = ft_volumerealign(cfg, t1, mri_acpc);
    % mri_t12acpc       = ft_volumerealign(cfg, t1, mri_acpc_brain);
    
    % save transform
    transform_t12ctf    = mri_t12acpc.transform; 
    save(fullfile(mripath, ...
                  sprintf('%s_transform_t12ctf_n',subj_match.pid{ss})), ...
                  'transform_t12ctf');

    % save registered volume
    cfg             = [];
    cfg.filename    = fullfile(mripath,sprintf('%s_t12ctf_n', subj_match.pid{ss}));
    cfg.filetype    = 'nifti';
    cfg.parameter   = 'anatomy';
    ft_volumewrite(cfg, mri_t12acpc);
    
    % freesurfer -> acpc -> ctf
    transform_surf2mri = mri.transform/mri_acpc.transform * ...
                         mri_t12acpc.transform/mri_t12acpc.transformorig;
   
    %%% LOAD FREESURFER SURFACE--------------------------------------------
    filename            = fullfile(datapath,[subj_match.pid{ss},'.L.midthickness.4k_fs_LR.surf.gii']);
    sourcemodel_4korig  = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});

    % transform to ctf space
    sourcemodel_4k      = ft_transform_geometry(transform_surf2mri, sourcemodel_4korig);
    sourcemodel_4k.inside = sourcemodel_4k.atlasroi>0;
    sourcemodel_4k      = rmfield(sourcemodel_4k, 'atlasroi');
    
    %%% CREATE HEAD MODEL---------------------------------------------------

    % load template_grid (can create your own template grid)
    % NOTE: the path to the template file is user-specific
    ftpath      = which('ft_defaults.m');
    ftpath      = ftpath(1:end-14);
    a           = load(fullfile(ftpath, '/template/sourcemodel/standard_sourcemodel3d10mm'));
    template_grid = a.sourcemodel;

    % segment (brain mask) the anatomical CTF MRI for singleshell model
    cfg                 = [];
    cfg.output          = {'brain','skullstrip'};
    cfg.viewresult      = 'yes';
    cfg.spmversion      = 'spm12';
    seg                 = ft_volumesegment(cfg, mri);

    % construct the volume conductor model (i.e. head model) for each subject   
    % this is optional, and for the purpose of this tutorial only required for
    % plotting, later on
    cfg                 = [];
    cfg.method          = 'singleshell'; %config.beamforming.headmodel.method;
    cfg.spmversion      = 'spm12';
    headmodel           = ft_prepare_headmodel(cfg, seg);
    headmodel           = ft_convert_units(headmodel,'cm');
    
    % save head model and surface (4k)
    save(fullfile(ssSubjPath(ss),sprintf('/%s_sourcemodel_4k',subj_match.pid{ss})), 'sourcemodel_4k');
    save(fullfile(ssSubjPath(ss),sprintf('/%s_hdm',subj_match.pid{ss})), 'headmodel');

     %%% LOAD MEG DATA ---------------------------------------------------------
    load([ssSubjPath(ss) '/ft_meg_fullyProcessed.mat'],'-mat','data');

    %% VISUALIZATION: Headmodel and Surfaces --------------------------------
    close all;
    
    figure;
    ft_plot_headmodel(headmodel, 'edgecolor', 'none'); alpha 0.5
    ft_plot_mesh(ft_convert_units(sourcemodel_4k, 'cm'),'vertexcolor',sourcemodel_4k.sulc);
    ft_plot_sens(data.grad);
    view([0 -90 0])
        
    % save fig
    saveas(gcf,fullfile(ssSubjPath(ss),'/surf_head.png'))
    savefig(fullfile(ssSubjPath(ss),'/surf_head.fig'))

end