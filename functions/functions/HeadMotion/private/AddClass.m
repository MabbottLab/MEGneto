% Add trial class.
function Class = AddClass(Class, Name, Trials)
  if ~isempty(Trials)
    nCl = numel(Class);
    Ix = nCl + 1; % Index for new trial class.
    for c = 1:nCl
      if strcmpi(Class(c).Name, Name)
        % Class already exists.
        Ix = c;
        break
      end
    end
    if Ix > nCl
      % Create trial class.
      Class(Ix).Name = Name;
      Class(Ix).Id = max([Class(:).Id]) + 1;
      %     nCl = nCl + 1;
    end
    Class(Ix).Trials = unique([Class(Ix).Trials; ...
      Trials(:)]);
    Class(Ix).Count = numel(Class(Ix).Trials);
  end
end
