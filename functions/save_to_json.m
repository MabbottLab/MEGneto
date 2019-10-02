function save_to_json(data, filepath, overwrite)
%SAVE_TO_JSON safer, more limited savejson wrapper that
%turns tables to structs for saving for Megne2
% save_to_json(data, filepath, overwrite)
if nargin < 3
    overwrite = false;
end
if isa(data,'table')
    data = table2struct(data);
end
if overwrite == false && exist(filepath, 'file')
    warning([filepath ' already exists and overwrite disabled. ' ...
        'Use existing config to modify values if required.']);
else
    savejson('',data,filepath);
end