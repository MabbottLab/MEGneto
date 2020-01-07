function path_check(paths, exclude)

% PATH_CHECK verify that all paths and folders were initialized correctly.
%
% INPUTS:
%   paths               =   struct returned by PATH_GENERATION 
%   exclude             =   cell-array of excluded fields
%
% RETURNS:
%   displays warning if file or folder does not exist.  
%
% See also: MEGNE2SETUP, PATH_GENERATION

% Last updated by: Julie Tseng, 2020-01-07
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.

paths = table2struct(paths);
paths = rmfield(paths,exclude);
paths = struct2cell(paths);
for ii = 1:length(paths)
    bb = paths{ii};
    if ~exist(bb,'file')
        warning([bb ' does not exist'])
    end
end