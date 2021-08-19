% Generate list of good trials from Class structure.
function GoodTrials = GetGoodTrials(Class, nT)
  nCl = numel(Class);
  GoodTrials = (1:nT)';
  for c = 1:nCl
    if strncmpi(Class(c).Name, 'BAD', 3)
      GoodTrials = setdiff(GoodTrials, Class(c).Trials);
    end
  end
end