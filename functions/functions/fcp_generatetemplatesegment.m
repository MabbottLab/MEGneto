function [ template_seg ] = fcp_generatetemplatesegment( )
% FCP_GENERATETEMPLATESEGMENT check if already-segmented template exists. if so, use it. otherwise, generate!

% determine path to template MRI
fullPath = which('ft_preprocessing.m');
if ~isempty(fullPath)
    [pathstr,~,~] = fileparts(fullPath);
    info.template_path = [pathstr,'/template/anatomy/single_subj_T1.nii'];
    info.template_coordsys = 'spm';        % set to [] to use the default coordinate system
else
    disp('fieldtrip cound not be found. Make sure it is in your path');
end

[~,file,~] = fileparts(info.template_path);
[pathstr,~,~] = fileparts(which('fcp_generatetemplatesegment'));

% check if already-segmented template exists
if ~exist([pathstr '/' file '.mat'], 'file')
    
    % Load MNI template brain
    template = ft_read_mri(info.template_path);
    template.coordsys = info.template_coordsys;
    
    % Segment the template brain
    cfg = [];
    cfg.output = {'brain', 'skull', 'scalp'};
    template_seg = ft_volumesegment(cfg,template);
    
    save([pathstr '/' file '.mat'], 'template_seg', '-v7.3');
else
    load([pathstr '/' file '.mat'], 'template_seg');
end



end

