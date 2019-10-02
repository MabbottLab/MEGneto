function path_check(paths, exclude)
%PATH_CHECK checks whether all the strings within the provided struct are
%valid existing files / directories, with fields specified within the cell
%array "excluded" not included
paths = table2struct(paths);
paths = rmfield(paths,exclude);
paths = struct2cell(paths);
for ii = 1:length(paths)
    bb = paths{ii};
    if ~exist(bb,'file')
        warning([bb ' does not exist'])
    end
end