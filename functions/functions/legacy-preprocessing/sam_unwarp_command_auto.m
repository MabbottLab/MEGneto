function coords = sam_unwarp_command_auto(sn3d, svl, coordsFile, talUnits)


% 
%   
%  
%   coords = sam_unwarp_command_auto(sn3d, svl, coordsFile, talUnits)
%  
%   [sn3d] is matrix of sn3d.mat files
%   [svl] is matrix of an unwarped a*img file for each subject
%   [coordsFile] is mx3 matrix of mni or tlc coords. m=how many points.
%   [talUnits]=1 for tlrc, 0 for mni
%   [coords]=tlrc or mni output coords

[m n] = size(svl);

for ii = 1:m
    coords(ii) = sam_unwarp_command(sn3d(ii,:), svl(ii,:), coordsFile, talUnits);
end
%
% something has changed with matlab and definding the field seems to tack
% another on top, below is what used to work.
%    coords(ii).subject = sam_unwarp_command(sn3d(ii,:), svl(ii,:), coordsFile, talUnits);