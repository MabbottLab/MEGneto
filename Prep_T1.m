function Prep_T1(config, pid, mri_path, visitnum)

% setup: load things
if exist('visitnum', 'var')
    this_output = fullfile(config.meta.project_path, ...
                           config.meta.analysis_name, ...
                           pid, ...
                           sprintf('ses-%.2d', visitnum));
else
    this_output = fullfile(config.meta.project_path, ...
                       config.meta.analysis_name, ...
                       pid);
end
load([this_output '/out_struct.mat']);


%% Load MRI

% make sure it's a char path
mri = ft_read_mri(mri_path);

% realign
cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'ctf';
mri = ft_volumerealign(cfg, mri);
% mark lpa, rpa, npa, and +z point

% save the result
save([this_output '/Prep_T1_aligned.mat'], 'mri');

end