% [MeanContrast, SigContrast, pValue] = PermutationTest( ...
%    Data, Data2, nPerms, pValue, Paired, OneSided, TMaxDims, ...
%    MultipleComparisons)
%
% One- or two-sample, paired or unpaired random permutation test of the
% null hypothesis that the samples in arrays Data and Data2 along the last
% dimension come from the same distribution(s). Samples along dimensions
% other than the last are either treated as independent tests (with
% possibility of multiple comparison correction) or an "global" test
% ("tmax" of Blair and Karniski, 1993) is performed.  A combination of
% independent and global tests can be performed on different dimensions of
% the data arrays.
%
%   *Inputs*
% Data, Data2: Data arrays representing 2 groups of subjects. The data must
% be arranged such that the last dimension represents the subjects, i.e.
% the individuals that will be permuted between groups. E.g. for 3-d voxel
% images, Data(:, :, :, 2) would be the data for the 2nd subject. If Data2
% is not provided, it is taken to be 0, i.e. a one-sample test is performed
% on Data.
%
% All other input arguments are optional, default values in square
% brackets.
%
% nPerms [2048]: Number of permutations to test. pValue [0.05]: The desired
% threshold of significance for the test.
%
% Paired [true if same number of subjects, false otherwise]: Whether to
% perform a paired or unpaired test.  Paired requires the same number of
% subjects in both groups (logically the same ones).
%
% OneSided [false]: Whether to apply the requested pValue to the positive
% side (tail) only, investigating only the contrast Data>Data2, or to both
% sides.  When asking for a two-sided test (OneSided=false), the
% distribution of absolute values is used, which is equivalent to the
% Bonferroni correction, dividing the pValue equally between the 2 sides.
%
% TMaxDims [[] (empty)]: Dimensions of the data arrays to be "grouped" in
% order to perform a tmax test.  This means we assume data points along
% these dimensions belong to the same distribution.  If this is not the
% case, strong effects will dominate and smaller effects that would have
% been significant if tested independently may not survive.
%
% MultipleComparisons [1 (no correction)]: Number of comparisons to take
% into account when applying the Sidak correction (assumes independent
% tests) such that the familywise error rate will be given by the provided
% pValue. Because in many cases the number of comparisons is not simply
% equal to the number of data points, e.g. for time samples it depends on
% the filtering bandwidth and for voxels it is at most equal to the number
% of sensors, this value is left to the user to provide. The default value
% is 1, i.e. no correction.
%
% Important.  While both previous options can be used simultaneously,
% MultipleComparisons is meant to be used for dimensions that are not
% included in TMaxDims, as the tmax test does not require multiple
% comparison correction.
%
%   *Outputs*
% MeanContrast: Averaged contrasted data between the two groups.
%
% SigContrast: MeanContrast where the pValues are smaller than the provided
% pValue threshold, i.e. where the contrast is significant. 
% 
% pValues: Array of false positive probability values (p) for each of the
% data points of the MeanContrast.
%
% Dependencies: d2b, b2d
%
% 2014-05-12
% Marc Lalancette, The Hospital for Sick Children, Toronto, Canada.
% Inspired by code from Travis Mills.

function [MeanContrast, SigContrast, pValues] = PermutationTest( ...
    Data, Data2, nPerms, pValue, Paired, OneSided, ...
    TMaxDims, MultipleComparisons)
  
  % Ensure different pseudo-randomness each run.
  RandomSettings = rng;
  if RandomSettings.Seed == 0
    rng('shuffle');
  end

  % -----------------------------------------------------------------------
  % Code testing stuff.
  
  % Test subfunctions.
  %   n = 55;
  %   k = round(n/2);
  %   Count = 50;
  %   Verify = true;
  %   [Sequences, Count] = RandomCombinations(n, k, Count, Verify)
  %   [Sequences, Count] = RandomSequences(n, Count, Verify);
  %   sort(sum(Sequences * 2.^(0:size(Sequences, 2)-1)', 2))
  %   sum(Sequences, 2)/n
  %   keyboard
  %   return
  
  %   nO = 10;
  %   nI = 1;
  %   nS = 4;
  %   Data = rand(nO, nI, floor(nS/2));
  %   Data2 = rand(nO, nI, ceil(nS/2));
  
  %   Data = [-0.0836079479005819	1.64905393375920	1.91558100652028	2.71956223234769	0.856750862847083	-0.496310537137662	0.841230689163967	1.11927882127676	0.520849595255211	-0.262047436272037	-0.529848900930624];
  %   Data2 = [0.755749469505206	0.436374974711638	0.752830391856942	2.00064904229199	-1.45872550554240	-1.07610802246434	0.0561047887419034	1.30120906140266	-0.272329774491732	-0.0655473066489307	-0.852122804066966];
  %   nPerms = 50000;
  %   pValue = 0.05;
  %   Paired = false;
  %   OneSided = false;
  %   TMaxDims = [];
  %   MultipleComparisons = 1;
  %   nargout = 3;
  
  
  % -----------------------------------------------------------------------
  % Get and test input arguments.  And organize data.
  
  if ~exist('OneSided', 'var') || isempty(OneSided)
    OneSided = true;
  end
  
  if ~exist('TMaxDims', 'var')
    TMaxDims = [];
  end
  
  if ~exist('MultipleComparisons', 'var') || isempty(MultipleComparisons)
    MultipleComparisons = 1;
  end
  
  if ~exist('pValue', 'var') || isempty(pValue)
    pValue = 0.05;
  end
  
  if ~exist('nPerms', 'var') || isempty(nPerms)
    nPerms = 2048;
  end
  
  % Everything is at least a 2-dim array to Matlab, and additional trailing
  % singleton dimensions are ignored. If there is 1 subject in one group,
  % one dimension will appear missing.
  nDataDims = max(ndims(Data), ndims(Data2));
  DataSize = size(Data);
  if ~exist('Data2', 'var') || isempty(Data2)
    Data2 = zeros(DataSize);
    if exist('Paired', 'var') && ~Paired
      fprintf('One-sample paired and unpaired tests are equivalent, using faster paired test.');
    end
    Paired = true;
  end
  DataSize2 = size(Data2);
  % Reject second dimension if column vector. This will only happen if
  % testing scalars.
  if nDataDims == 2 && DataSize(2) == 1 && DataSize2(2) == 1
    nDataDims = nDataDims - 1;
    DataSize(2) = [];
    DataSize2(2) = [];
  end
  % Check data matches between groups.
  if any(DataSize(1:nDataDims-1) ~= DataSize2(1:nDataDims-1))
    error('Data array size mismatch between groups.')
  end
  % Verify that tmax dimensions correspond to data and not subjects.
  if any(TMaxDims >= nDataDims)
    error('Requested TMaxDims exceed data dimensions, assuming the last dim corresponds to subjects');
  end
  
  nDs1 = DataSize(nDataDims);
  nDs2 = DataSize2(nDataDims);
  nDatasets = nDs1 + nDs2;
  if nDatasets < 3
    error('You need at least 2 subjects in one group to permute.');
  end
  Asymmetrical = nDs1 ~= nDs2;
  
  if ~exist('Paired', 'var') || isempty(Paired)
    if Asymmetrical
      Paired = false;
      fprintf('Performing unpaired permutation test.\n');
    else
      Paired = true;
      fprintf('Performing paired permutation test.\n');
    end
  end
  
  if Paired && Asymmetrical
    error('Paired permutation test requires the same number of subjects in both groups, i.e. data arrays should have the same size.')
  end
  
  % Combine the data in one array.
  Data = cat(nDataDims, Data, Data2);
  clear Data2 DataSize2;
  DataSize(nDataDims) = nDatasets;
  
  
  % -----------------------------------------------------------------------
  % pValue corrections and checks.
  % Sidak "inverse correction" for multiple comparisons.
  RequestedpValue = pValue;
  pValue = MultiTestCorrection(RequestedpValue, MultipleComparisons, 1, true);
  %   pValue = 1 - (1 - pValue)^(1/MultipleComparisons);
  
  % Check if requested significance level is achievable.
  if Paired
    MaxPerms = 2^nDs1;
  else
    MaxPerms = nchoosek(nDatasets, nDs1);
  end
  if nPerms >= MaxPerms
    nPermsTemp = MaxPerms;
  else
    nPermsTemp = nPerms;
  end
  if ~OneSided && ~Asymmetrical
    % The Bonferroni correction for two-sided test is taken into account by
    % taking the absolute value below.  There is no need to change the
    % pValue here.  But we do end up effectively with half as many points
    % in the permutation distribution since each point will appear twice.
    nPermsTemp = ceil(nPermsTemp/2);
  end
  if nPermsTemp < 1/pValue
    % Check if maximum number of permutations for number of subjects.
    if nPerms >= MaxPerms
      error(['The requested (possibly corrected) significance level (pValue=%1.2g) cannot be achieved with the provided number of subjects. ', ...
        'Either reduce the number of independent data points and correspondingly the number of tests (MultipleComparisons), ', ...
        'or consider doing a TMax or uncorrected test.'], RequestedpValue);
    else
      error(['The requested (possibly corrected) significance level (pValue=%1.2g) requires at least %d permutations '...
        '(which is possible with the provided number of subjects).'], RequestedpValue, ceil(1/pValue));
    end
  else
    CorrectedPrecision = MultiTestCorrection(1/nPermsTemp, MultipleComparisons, 1, false) / RequestedpValue;
    if CorrectedPrecision > 0.05
      if nPermsTemp >= MaxPerms
        warning(['The number of permutations cannot be increased and would give a poor precision (%1.0f%%) ', ...
          'in terms of the (possibly corrected) significance level (pValue=%1.2g). ', ...
          'It is recommended to either reduce the number of independent data points and correspondingly the number of tests (MultipleComparisons), ', ...
          'or consider doing a TMax or uncorrected test.'], ...
          CorrectedPrecision * 100, RequestedpValue);
      else
        warning(['The requested number of permutations would give a poor precision (%1.0f%%) ', ...
          'in terms of the (possibly corrected) significance level (pValue=%1.2g). ', ...
          'It is recommended to increase the number of permutations.'], ...
          CorrectedPrecision * 100, RequestedpValue);
      end
    end
  end
  
  
  % -----------------------------------------------------------------------
  % Reshape data to simplify and speed up testing.  First dimension will be
  % "grouped" data from TMaxDims.  Second dimension will be independent
  % data dimensions.  Third: subjects.  Use singleton dimensions if no
  % tmax test or no independent data points.
  %   TMaxDims = sort(TMaxDims(DataSize(TMaxDims) > 1));
  %   Singletons = setdiff(1:nDataDims-1, TMaxDims);
  %   NonTMaxDims = Singletons(DataSize(Singletons) > 1);
  %   Singletons = setdiff(Singletons, NonTMaxDims);
  %   DimPermutation = [TMaxDims, NonTMaxDims, nDataDims, Singletons];
  
  Singletons = [find(DataSize == 1), nDataDims + (1:2)];
  DataSize(end+(1:2)) = 1;
  if isempty(TMaxDims)
    % Use first singleton dimension.
    TMaxDims = Singletons(1);
  else
    % Ensure dimensions are sorted so we can get the correct shape again.
    TMaxDims = sort(TMaxDims);
  end
  NonTMaxDims = setdiff(1:nDataDims-1, TMaxDims);
  if isempty(NonTMaxDims)
    % Use first singleton dimension, other than those in TMaxDims.
    Singletons = setdiff(Singletons, TMaxDims);
    NonTMaxDims = Singletons(1);
  end
  DimPermutation = [TMaxDims, NonTMaxDims, nDataDims];
  [~, InverseDimPerm] = sort(DimPermutation);
  nTMaxElements = prod(DataSize(TMaxDims)); % prod([]) = 1
  nIndepElements = prod(DataSize(NonTMaxDims));
  Data = reshape(permute(Data, DimPermutation), ...
    nTMaxElements, nIndepElements, nDatasets, []);
  
  % We now use bsxfun which is faster than indexing with ones like repmat.
  %   OnesVectors = {ones(1, nTMaxElements), ones(1, nIndepElements)}; 
  
  
  % -----------------------------------------------------------------------
  % Main section.
  
  % Compute the real contrast first.
  % This formula will do proper averages before contrast even if groups
  % have different numbers in unpaired condition.
  %   Coefficients(OnesScalar{:}, :) = [ones(1, nDs1)./nDs1, -ones(1, nDs2)./nDs2];
  %   MeanContrast = sum(Data .* Coefficients(OnesVectors{:}, :), nDataDims);
  Coefficients(1, 1, :) = [ones(1, nDs1)./nDs1, -ones(1, nDs2)./nDs2];
  MeanContrast = sum(bsxfun(@times, Data, Coefficients), 3);
  
  % Get rid of values that are due to machine precision errors (make them zero).
  % Total error in contrast values.
  Precision = sum(eps(Data), 3);
  MeanContrast(abs(MeanContrast) <= Precision) = 0;
  
  % Stop here if not asking for anything else.
  if nargout < 2
    % Give back original data shape.
    MeanContrast = permute( reshape(MeanContrast, DataSize(DimPermutation(1:end-1))), InverseDimPerm );
    return;
  end
  
  % Generate the permutations as logical arrays indicating which group each
  % subject belongs to.
  if Paired
    % Permutation distribution is symmetrical so we can do only half as
    % many permutations.  For the one-sided test, we get 2 different points
    % from each permutation (-min and max), while for the two-sided test,
    % we use max(abs()), effectively the Bonferroni correction, which would
    % give twice each result over all permutations.
    %
    % So instead, we set the first subject in the first group and permute
    % the rest only.  The "negative" side then corresponds to all
    % permutations where the first subject would be in the second group.
    [Combinations, nPerms] = RandomSequences(nDs1 - 1, nPerms/2);
    if OneSided
      % We get 2 points for each permutation for the distribution.
      nDistrib = 2 * nPerms;
    else
      % We only get 1 point per permutation.
      nDistrib = nPerms;
    end
    Combinations = [true(nPerms, 1), Combinations, ...
      false(nPerms, 1), ~Combinations];
  else % unpaired.
    if ~Asymmetrical
      % Distribution is symmetrical as for the paired test; only use half.
      [Combinations, nPerms] = RandomCombinations(nDatasets - 1, nDs1 - 1, nPerms/2);
      Combinations = [true(nPerms, 1), Combinations];
      if OneSided
        % We get 2 points for each permutation for the distribution.
        nDistrib = 2 * nPerms;
      else
        % We only get 1 point per permutation.
        nDistrib = nPerms;
      end
      
    else % Different number of subjects in the two groups.
      % "Opposite contrasts" are no longer part of the distribution so we
      % must use all nPerms permutations and look only at max for 1-sided
      % test.
      [Combinations, nPerms] = RandomCombinations(nDatasets, nDs1, nPerms);
      nDistrib = nPerms;
    end
  end
  
  % Get the distributions of contrast values across permutation of subjects
  % between the two groups, for each independent (i.e. non-grouped for the
  % tmax test) data point.
  Distrib(nIndepElements, nDistrib) = 0;
  for p = 1:nPerms
    Coefficients(Combinations(p, :)) = 1/nDs1;
    Coefficients(~Combinations(p, :)) = -1/nDs2;
    Contrast = sum(bsxfun(@times, Data, Coefficients), 3);
    % Zero values smaller than precision.
    Contrast(abs(Contrast) <= Precision) = 0;
    if nTMaxElements > 1
      % Find max (over tmax group) distribution for each data point.
      % (Best always to give dimension explicitly to max and min in case of
      % row vector, though we just tested for this.)
      if OneSided
        Distrib(:, p) = max(Contrast, [], 1); % max over first dim.
        if ~Asymmetrical
          Distrib(:, p + nPerms) = -min(Contrast, [], 1);
        end
      else % 2-sided, use absolute value.
        Distrib(:, p) = max(abs(Contrast), [], 1); % max over first dim.
      end
    else % No tmax test, no need for max.
      % Find distribution for each data point.
      if OneSided
        Distrib(:, p) = Contrast;
        if ~Asymmetrical
          Distrib(:, p + nPerms) = -Contrast;
        end
      else
        Distrib(:, p) = abs(Contrast);
      end
    end
  end
  % Sort the distributions.
  Distrib = sort(Distrib, 2);
  
  
  % Find significant values in mean contrast.
  % Index along distribution where the significance threshold should be.
  ThresholdIndex = floor((1 - pValue) * nDistrib) + 1;
  Thresholds(1, :) = Distrib(:, ThresholdIndex);
  % Check for duplicate value of the chosen threshold(s) preceding it in
  % the distribution, which would require finding the next larger value(s)
  % to use as threshold(s).
  if ThresholdIndex > 1 % (1 is possible if pValue = 1 or almost.)
    if any( Thresholds(:) == Distrib(:, ThresholdIndex - 1))
      NeedFixing = find(Thresholds(:) == Distrib(:, ThresholdIndex - 1));
      for i = NeedFixing(:)' % Older Matlab requires row vector.
        t = ThresholdIndex + 1;
        while t <= nDistrib && Thresholds(i) == Distrib(i, t)
          t = t + 1;
        end
        if t > nDistrib
          % No larger value found.  Nothing will be above threshold.
          Thresholds(i) = Inf;
        else
          % Found new threshold.
          Thresholds(i) = Distrib(i, t);
        end
      end
    end
  end
  
  if OneSided
    SigIndices = find( bsxfun(@gt, MeanContrast, Thresholds) ); % MeanContrast > Thresholds(OnesVectors{1}, :)
  else
    SigIndices = find( bsxfun(@gt, abs(MeanContrast), Thresholds) ); % abs(MeanContrast) > Thresholds(OnesVectors{1}, :)
  end
  
  % Create array with only significant values.
  SigContrast = NaN(nTMaxElements, nIndepElements);
  SigContrast(SigIndices) = MeanContrast(SigIndices);
  
  % Calculate pValues for each data point of the mean contrast.
  if nargout > 2
    pValues(nTMaxElements, nIndepElements) = 0;
    for i = 1:nIndepElements
      for o = 1:nTMaxElements
        if OneSided
          Index = find(Distrib(i, :) >= MeanContrast(o, i), 1, 'first');
        else
          Index = find(Distrib(i, :) >= abs(MeanContrast(o, i)), 1, 'first');
        end
        if isempty(Index)
          % This can happen when not doing every possible permutation.
          % Would give a pValue of 0, so return instead 0 + uncertainty.
          Index = nDistrib;
        end
        pValues(o, i) = 1 - (Index - 1)/nDistrib;
      end
    end
    % Sidak correction for multiple comparisons.
    pValues = 1 - (1 - pValues).^(MultipleComparisons);
    % Give back original data shape.
    pValues = permute( reshape(pValues, DataSize(DimPermutation(1:end-1))), InverseDimPerm );
  end
  
  % Give back original data shape.
  MeanContrast = permute( reshape(MeanContrast, DataSize(DimPermutation(1:end-1))), InverseDimPerm );
  SigContrast = permute( reshape(SigContrast, DataSize(DimPermutation(1:end-1))), InverseDimPerm );
  
  % For testing.
  %   figure;
  %   hist(Distrib, -2:0.1:2);
  %   pValues
end


%------------------------------------------------------------------------
% This is used for paired permutation test.
function [Sequences, Count] = RandomSequences(n, Count, Verify)
  % 'Sequences' here means combinations for any number of choices (k) among n.
  % Sequences (size = [Count, n]) is a logical array where 'chosen' elements are 'true'.
  
  % As n increases, there is often little chance of duplicates, so we can
  % probably omit verification for speed increase if a few duplicates are ok.
  if ~exist('Verify', 'var') || isempty(Verify)
    Verify = true;
  end
  
  Verbose = false;
  if Count >= 2^n
    Count = 2^n;
    Exact = true;
    if Verbose
      fprintf('Exact test with %d permutations.\n', Count);
    end
  else
    Exact = false;
    if Verbose
      fprintf('Random test with %d permutations.\n', Count);
    end
  end
  
  if eps(2^n) < 1 % n < 52 % eps(2^52) = 1 and is therefore no longer accurate as an integer.
    % Use first method: binary sequence.
    % This is the fastest method, even more than filling directly (like
    % in combn_mod).
    %     Sequences = false(n, Count);
    if Exact
      % All possible combinations for any number of choices among n.
      Picks = 0:2^n-1; % = zeros(100, 'uint64');
    else
      if verLessThan('matlab', '7.13.0')
        % Older Matlab doesn't have randperm(n,k), and randperm(2^n)
        % required too much memory.
        Picks = my_randperm(2^n, Count) - 1; % Random distinct integers, not sorted.
      else
        Picks = randperm(2^n, Count) - 1; % Random distinct integers, not sorted.
      end
    end
    Sequences = d2b(Picks', n); % Logical array of binary digits, least significant first.
    %     for bit = 1:n
    %       Sequences(bit, :) = logical(bitget(Picks, bit)); % Least significant bit first.
    %     end
    
  else
    % Method two: random construction, need to check for duplicates.
    if Exact || Count > 2^(n-1) % Could verify memory demand for large Count too, prevent out of memory errors.
      % Random picks might generate many duplicates, but also way too many.
      error('Asked for too many permutations: %f', Count);
    end
    Sequences = rand(Count, n) < 0.5;
    % Verify and replace duplicate sequences, but duplicates are almost
    % impossible with reasonable numbers.
    if Verify
      s = 2;
      while s < Count
        if ~all( any( ...
            bsxfun(@xor, Sequences(1:s, :), Sequences(s+1, :)) ... % xor(Sequences(1:s, :), ones(s, 1)*Sequences(s+1, :))
            , 2) ) % if not all different, any not the same between combinations
          Sequences(s+1, :) = rand(1, n) < 0.5;
        else
          s = s + 1;
        end
      end
    end
    
  end % method choice
end % function RandomSequences


%------------------------------------------------------------------------
% This is used for unpaired permutation test.
function [Combinations, Count] = RandomCombinations(n, k, Count, Verify)
  % Choose k elements among n.
  % Combinations (size = [Count, n]) is a logical array where 'chosen' elements are 'true'.
  
  % As n increases, there is often little chance of duplicates, so we can
  % probably omit verification for speed increase if a few duplicates are ok.
  if ~exist('Verify', 'var') || isempty(Verify)
    Verify = true;
  end
  % Avoid warning from nchoosek precision limits if numbers are too big anyway.
  if eps(2^n) < 1
    MaxCount = nchoosek(n, k);
  else
    MaxCount = 1/eps('double');
  end;
  if Count >= MaxCount
    Count = MaxCount;
    if eps(Count) >= 1
      error('Asked for too many permutations: %f', Count);
    end
    Exact = true;
    fprintf('Exact test with %d permutations.\n', Count);
  else
    fprintf('Random test with %d permutations.\n', Count);
    if Count >= MaxCount/2
      if eps(MaxCount) >= 1
        error('Asked for too many permutations: %f', Count);
      end
      % Do exact then randomly pick the rows to keep.  Faster than trying to
      % randomly pick most of the possible different combinations (many
      % duplicates to correct).
      if verLessThan('matlab', '7.13.0')
        % Older Matlab doesn't have randperm(n,k), and randperm required
        % much memory.
        CombIndices = my_randperm(MaxCount, Count); % The combinations that we will actually keep.
      else
        CombIndices = randperm(MaxCount, Count); % The combinations that we will actually keep.
      end
      Count = MaxCount;
      Exact = true;
    else
      if eps(Count) >= 1
        error('Asked for too many permutations: %f', Count);
      end
      Exact = false;
    end
  end
  
  Combinations(Count, n) = false;
  if Exact
    Picks = nchoosek(1:n, k);
    for c = 1:Count % combination
      Combinations(c, Picks(c, :)) = true;
    end
    if exist('CombIndices', 'var')
      Combinations = Combinations(CombIndices, :);
      Count = length(CombIndices);
    end
  else
    
    % Haven't thought of a fast non-looping method yet, and that wouldn't
    % require verification.
    nChosen = zeros(Count, 1);
    for e = 1:n % element
      Combinations(:, e) = rand(Count, 1) < ((k - nChosen) ./ (n - e + 1));
      nChosen = nChosen + Combinations(:, e);
    end
    
    % Verify and replace duplicate combinations.
    if Verify
      c = 1;
      if eps(2^n) < 1
        % Use much faster binary method.
        CombLabels = b2d(Combinations);
        while c < Count
          if ismember(CombLabels(c+1), CombLabels(1:c)) % Same as any previous combination.
            nChosen = 0;
            for e = 1:n % element
              Combinations(c+1, e) = rand(1) < ((k - nChosen) ./ (n - e + 1));
              nChosen = nChosen + Combinations(c+1, e);
            end
            CombLabels(c+1) = b2d(Combinations(c+1, :));
          else
            c = c + 1;
          end
        end
      else
        % The if expression here is rather slow. About 25% is due to
        % replicating the combination.
        while c < Count
          % Maybe can use binary representation if eps(2^n) < 1
          if any( ~any( ...
            bsxfun(@xor, Combinations(1:c, :), Combinations(c+1, :)) ... % xor(Combinations(1:c, :), Combinations(ones(c, 1)*c+1, :) )
            , 2) ) % any previous combination witout any difference (same)
            nChosen = 0;
            for e = 1:n % element
              Combinations(c+1, e) = rand(1) < ((k - nChosen) ./ (n - e + 1));
              nChosen = nChosen + Combinations(c+1, e);
            end
          else
            c = c + 1;
          end
        end
      end
    end
    
  end % method choice
end % function RandomCombinations

%------------------------------------------------------------------------
% This is for older Matlab where randperm(n,k) doesn't exist.
function Perm = my_randperm(n, k)
  % Return k random distinct integers from the set 1:n.
  
  % If n is large and k is small, it takes too much memory to use the old
  % randperm(n).
  if k < n/2
    if verLessThan('matlab', '7.7.0')
      Perm = unique(ceil(rand(1, k) * n));
    else
      Perm = unique(randi(n, [1, k]));
    end
    % Replace duplicates that were removed by unique().
    while numel(Perm) < k
      Perm = unique([Perm, ceil(rand(1, k - numel(Perm)) * n)]);
    end
  else
    % k comparable to n, just truncate a full permutation.
    %     Perm = randperm(n);
    [unused, Perm] = sort(rand(1,n)); %#ok<ASGLU> % This is the old randperm(n).
    Perm(k+1:end) = [];
  end
end