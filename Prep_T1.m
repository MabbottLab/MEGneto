function Prep_T1(mri_path, out_path, pid, visitnum)

% Prep_T1 will load an MRI T1 into MATLAB and use a FieldTrip helper function 
% to make the T1 compatible for MEG beamforming (i.e., marking the 3 fiducial
% points (nasion, left, and right pre-auricular points). This only needs to be
% carried out once for each participant in the study, which is why the
% output should be saved to its own derivatives folder (e.g.,
% studies/derivatives/T1w_meg-beamforming), then referenced by the
% beamforming script. 
%
% INPUTS:
%   > mri_path:
%       string, full path to T1 file for this participant
%       NOTE: to use the subject-specific Glasser atlas parcellation later
%       on, you may want to use the post-FreeSurfer T1 for this step
%   > out_path:
%       derivatives folder where the output should be saved
%   > pid:
%       string, participant ID used to build the paths
%       to relevant data files and output folders
%   > visitnum:
%       int, visit number if longitudinal, optional arg)
%
% OUTPUTS:
%   > Prep_T1_aligned.mat:
%       *.mat file containing T1 MRI with fiducial locations marked, ready
%       for import into beamforming step
%
% See also: FT_VOLUMEREALIGN

% Last updated by: Julie Tseng, 2024-09-09
%   This file is part of MEGneto, see https://github.com/MabbottLab/MEGneto
%   for the documentation and details.

%% Load MRI

% make sure it's a char path
mri = ft_read_mri(char(mri_path));

% realign
cfg             = [];
cfg.method      = 'interactive';
cfg.coordsys    = 'ctf';
mri             = ft_volumerealign(cfg, mri);
% mark lpa, rpa, npa, and +z point

% make output directory if it doesn't already exist
full_outpath = fullfile(out_path, pid, sprintf("ses-%.2d", visitnum));
if ~exist(full_outpath, 'dir') 
    mkdir(full_outpath)
end

% save the result
save(fullfile(full_outpath, "Prep_T1_aligned.mat"), 'mri');

end