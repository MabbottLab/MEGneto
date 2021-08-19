function [ path ] = findPath( dirpath, pattern )
%FINDPATH Find a single file based on a given (regular expression?) pattern

if isunix
  [~,path] = unix(['find "',dirpath,'" -maxdepth 1 -iname "',pattern,'"']);
  % return the first result
  path = strsplit(path, '\n');
  if length(path) > 2
    warning('More than one path matched!\n%s', sprintf('%s\n',path{:}));
  end
  path = path{1};
else
  path = getfield(dir(fullfile(dirpath,pattern)),'name');
end


end

