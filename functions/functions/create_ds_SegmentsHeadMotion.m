function create_ds_SegmentsHeadMotion (ds_name, x, period, output)

% -------------------------------------------------------------------------
% Create new MarkerFile.mrk

% Check if markers are given as a vector
if ~isvector(x),
    error('Markers should be a vector.');
else 
    x = reshape(x,[],1);
end

% Check if DS file exists
if ~exist(ds_name),
    error('Source DS file does not exist.')
end

% make trl
cfg = [];
cfg.trl = [x, x+period-1, zeros(length(x),1)];
save(output, 'cfg');

% save([ds_name, '/ft_meg_grad_cfg.mat'], 'cfg');

%--------------------------------------------------------------------------
% Correct the initial head position
% HeadMotionTool('Dataset',ds_name,'CorrectInitial',true,'GUI',false);
                
return
