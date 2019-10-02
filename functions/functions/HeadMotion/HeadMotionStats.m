function [DistanceQ, MovementQ, TrialCounts, FitError] = ...
    HeadMotionStats(Datasets, MotionThreshold, SampleProportion, ...
    FitErrorThreshold, Verbose)
  % Calculate head motion statistics for a CTF MEG dataset.
  %
  % [Distance, Movement, TrialCounts, FitError] = ...
  %     HeadMotionStats(Datasets, MotionThreshold, SampleProportion, ...
  %     FitErrorThreshold, Verbose)
  %
  % Returns motion statistics to be used in analysis (e.g. for regressing
  % out potential motion effects from results, or compare different
  % populations).  All distances are in mm units.  This function should be
  % used after initial position correction and trial rejection with
  % HeadMotionTool.  This program uses the same computations as
  % HeadMotionTool, but does not modify datasets in any way.  It uses the
  % dataset's initial position and saved trial classification: only "good"
  % trials are used in calculating statistics.  The program will calculate
  % how many trials would be considered "bad" based on the provided
  % arguments (thresholds, etc.), but only trials already classified as
  % "bad" are not included in the statistics.  In other words, even trials
  % that would be rejected are included in e.g. the max distance, so it can
  % be greater than the provided threshold.
  %
  % Input arguments (Default values in brackets):
  %
  %  Datasets []: Cell array of strings: full path to dataset directories.
  %
  %  MotionThreshold [5]: Limit on movement distance from the initial
  % position (TrialsCount.BadDistance) or range of motion within a trial
  % (TrialsCount.BadMotion) for it to be considered "bad".  Depends on the
  % proportion of samples within the trial that are above this threshold
  % and the selected SampleProportion (see below).
  %
  %  SampleProportion [10]: Trials containing this proportion (in %) or more
  % of samples exceeding the thresholds will be considered "bad", i.e. the
  % quantile of the trial sample distribution used for thresholding.
  %
  %  FitErrorThreshold [10]: Trials where the coil localization fitting
  % error (in '%', greatest among the 3 coils) is above this threshold are
  % considered "bad" (TrialsCount.BadFit).
  %
  % Ouptuts: Each output variable is a structure and each field is a column
  % vector with one element for each input dataset.  The Mean field is the
  % trial average, Max is the maximum single-trial value.
  %
  %  DistanceQ: Mean, Max, Quantile.  Distance from reference among trials.
  %  This trial distance is the selected quantile of the sample distance
  %  distribution.  E.g. if SampleProportion = 10, the trial distance is
  %  the distance exceeded by 10% of samples within that trial. So one
  %  could say e.g. that 90% of trials have 90% of samples within the
  %  "Quantile" distance from the reference position, or that all trials
  %  have 90% of samples within the "Max" distance from the reference.
  %
  %  MovementQ: Mean, Max, Quantile, Total.  Range of motion within a trial
  %  (or for the whole dataset for Total).  Similarly to DistanceQ, this
  %  motion range is based on the selected proportion of samples.  (It is
  %  not calculated as precisely as distance, based on coil coordinates
  %  instead of a "rigid head".)
  %
  %  TrialCounts: BadDistance, BadMotion, BadFit, Rejected, Total.  Number
  %  of trials in each category.  "Rejected" indicates trials that were
  %  already classified as bad in the dataset.  The "bad" categories are
  %  based on the provided thresholds and only include trials that were not
  %  already rejected.  For a dataset processed with HeadMotionTool, one
  %  would expect zero "bad" trials, but possibly some already rejected.
  %
  %  FitError: Max.  Coil localization error, again at the selected
  %  proportion of samples within each trial.
  %
  %
  % Dependencies: This function and HeadMotionTool depend on a number of
  % files, some in a 'private' directory, and others, including the
  % unofficial CTF Matlab package from MISL, provided with SPM8 or
  % Fieldtrip, last updated 2012-04-16 (Linux specific bug fix).
  %
  % Marc Lalancette, The Hospital for Sick Children, Toronto, Canada.
  % 2014-03-25
  
  
  % Check some dependencies (not all, but should catch most cases).
  % See if dependencies should be added to the Matlab path.
  Dir = fileparts(mfilename('fullpath'));
  if ~exist('readCTFds.m', 'file')
    addpath([Dir, filesep, 'ctf']);
  end
  if ~exist('opencls.m', 'file')
    addpath(Dir);
  end
  if ~exist('readCTFds.m', 'file') || ~exist('opencls.m', 'file')
    error('Missing dependencies.  Make sure CTF Matlab files and other needed functions are in your Matlab path.');
  end
  
  if nargin < 1
    error('No datasets to process.');
  end
  if ~exist('MotionThreshold', 'var') || isempty(MotionThreshold)
    MotionThreshold = 5; % mm
  end
  if ~exist('SampleProportion', 'var') || isempty(SampleProportion)
    SampleProportion = 10; % percent
  end
  if ~exist('FitErrorThreshold', 'var') || isempty(FitErrorThreshold)
    FitErrorThreshold = 10; % percent
  end
  if ~exist('Verbose', 'var') || isempty(Verbose)
    if nargout == 0
      Verbose = true;
    else
      Verbose = false;
    end
  end
  
  handles.Phi = (1 + sqrt(5))/2;
  handles.Proportion = SampleProportion;
  handles.DistThresh = MotionThreshold;
  handles.FitThresh = FitErrorThreshold;
  
  if ischar(Datasets)
    Datasets = {Datasets};
  elseif ~iscell(Datasets)
    error('Unrecognized list of datasets, should be cell array.');
  end
  nD = numel(Datasets);
  
  Z = nan(nD, 1);
  DistanceQ = struct('Mean', Z, 'Max', Z, 'Quantile', Z);
  MovementQ = struct('Mean', Z, 'Max', Z, 'Quantile', Z, 'Total', Z);
  FitError = struct('Max', Z);
  TrialCounts = struct('BadDistance', Z, 'BadMotion', Z, 'BadFit', Z, ...
    'Rejected', Z, 'Total', Z);
  
  if Verbose
    clc
    fprintf('Units: mm\n');
    fprintf('Datasets:\n');
    for d = 1:nD
      fprintf(' %d: %s\n', d, Datasets{d});
    end
    fprintf('\n');
    fprintf('Dataset, Distance (mean, max, quantile), Movement (mean, max, quantile, total), Fit error (max), Trials (bad distance, bad motion, bad fit, already rejected, total)\n');
  end
  
  
  for d = 1:nD
    
    % Get data.
    handles.Dataset = Datasets{d};
    if ~isdir(handles.Dataset)
      warning('Dataset not found: %s\n', handles.Dataset);
      continue;
    end
    %   [handles.Fiducials, ...
    %     handles.Data, handles.FitData, handles.nS, handles.nT, ...
    %     handles.Class, handles.GoodTrials] = ...
    %     GetData(handles.Dataset);
    handles = GetData(handles);
    % Everything in dewar coordinates, in mm.
    % MedianLoc, Fids: [1, nC, nT]
    % Data: [nS, nC, nT]
    
    % Get rigid reference body, in head coordinates.
    handles = ReferenceBody(handles);
    
    
    % ----------------------------------------------------------------------
    % Get additional information on movement overall and per trial, possibly
    % for automated trial rejection, or to get values that can be used in
    % correlations or as weights in group averaging for example.
    nT = handles.nT;
    Reference = handles.Fiducials;
    
    % Calculate distance from chosen location (median or fids) across trials
    % and time.  Use geometric median location per trial for unique trial
    % distance. Note that this is not an actual position visited in the
    % trial, so it could be outside the range of MinDistance to MaxDistance.
    %     T_Location = GeoMedian(reshape(handles.Data, [handles.nS, 3, 3, nT]), 1e-3); % [1, 3, 3, nT]
    % Both RigidDistances and CorrectRigid accept [nS, 3, 3, nT] or [nS, 9, nT] shapes.
    
    % No need to correct for rigid shape here.
    %     T_LocationDist = squeeze(RigidDistances(T_Location, Reference)); % squeeze([1, 1, nT])
    
    % Distance and FitError quantile of samples per trials at the desired
    % proportion.
    T_P = max(1, min(handles.nS, round((1 - ...
      handles.Proportion/100 ) * handles.nS)));
    % Divide by two for range quantile, since we "cut" from both sides.
    T_P2 = max(1, min(handles.nS, round((1 - ...
      handles.Proportion/100 /2 ) * handles.nS)));
    % Only use data from good trials for overall measures.
    nTG = length(handles.GoodTrials);
    if nTG > 0 % In case of all trials bad.
      
      % min and max distances and desired quantile.
      Data = CorrectRigid(handles.Data, handles.Rigid);
      Data = bsxfun(@minus, Data, Reference);
      Distances = RigidDistances(handles.Data, Reference); % [nS, 1, nT]
      
      T_Distance_Q = sort(Distances, 1);
      %       T_MaxDistance = squeeze(T_Distance_Q(handles.nS, :, :)); % squeeze([1, 1, nT])
      %     T_MinDistance = squeeze(T_Distance_Q(1, :, :));
      T_Distance_Q = squeeze(T_Distance_Q(T_P, :, :)); % nT
      
      DistanceQ.Mean(d) = mean(T_Distance_Q(handles.GoodTrials));
      DistanceQ.Max(d) = max(T_Distance_Q(handles.GoodTrials));
      TG_P = max(1, min(nTG, round((1 - handles.Proportion/100 ) * nTG)));
      TDQ_Ordered = sort(T_Distance_Q(handles.GoodTrials));
      DistanceQ.Quantile(d) = TDQ_Ordered(TG_P);
      
      % Upper bound for range of movement.  Max out of the 3 coils of sum of
      % max range in each direction.
      Data = sort(Data, 1);
      %       T_Range = squeeze(MaxCoilNorms(Data(end, :, :) - Data(1, :, :))); % squeeze([1, 1, nT])
      T_Range_Q = squeeze(MaxCoilNorms(Data(T_P2, :, :) - Data(end+1-T_P2, :, :))); % squeeze([1, 1, nT])
      % Real range can't exceed twice the max distance so enforce that
      % bound.  To be a bit less conservative, use distance quantile for
      % range quantile.
      %       T_Range = min(T_Range, 2 * T_MaxDistance);
      T_Range_Q = min(T_Range_Q, 2 * T_Distance_Q);
      MovementQ.Mean(d) = mean(T_Range_Q(handles.GoodTrials));
      MovementQ.Max(d) = max(T_Range_Q(handles.GoodTrials));
      TRQ_Ordered = sort(T_Range_Q(handles.GoodTrials));
      MovementQ.Quantile(d) = TRQ_Ordered(TG_P);
      
      % Quantile of the whole sample*trial distribution.
      P = max(1, min(handles.nS*nTG, round((1 - ...
        handles.Proportion/100 ) * handles.nS*nTG)));
      P2 = max(1, min(handles.nS*nTG, round((1 - ...
        handles.Proportion/100 / 2) * handles.nS*nTG)));
      DistancesG = Distances(:, :, handles.GoodTrials);
      Distance_Q = sort(DistancesG(:));
      Distance_Q = Distance_Q(P);
      Range = reshape(permute( Data(:, :, handles.GoodTrials), ...
        [1, 3, 2]), [handles.nS * nTG, 9]);
      Range = sort(Range, 1);
      Range_Q = MaxCoilNorms(Range(P2, :) - Range(end+1-P2, :)); % squeeze([1, 1])
      %       Range = MaxCoilNorms( max(Range, [], 1) - min(Range, [], 1) );
      Range_Q = min(Range_Q, 2 * Distance_Q);
      MovementQ.Total(d) = Range_Q;
      
      
      % ----------------------------------------------------------------------
      % Average coil position fitting error, convert to percent.  Always use
      % max of any coil direction.  Only changes if changing thresholds.
      FitData = max(100 * handles.FitData, [], 2); % [nS, 1, nT]
      %     T_Fit = squeeze(mean(FitData, 1)); % squeeze([1, 1, nT])
      %     Fit = mean(T_Fit(handles.GoodTrials));
      
      T_Fit_Q = sort(FitData, 1);
      %     T_MaxFit = squeeze(T_Fit_Q(end, :, :)); % nT
      %     T_MinFit = squeeze(T_Fit_Q(1, :, :)); % nT
      T_Fit_Q = squeeze(T_Fit_Q(T_P, :, :)); % nT
      
      %       FitData = FitData(:, :, handles.GoodTrials);
      %   Fit_Q = sort(FitData(:));
      FitError.Max(d) = max(T_Fit_Q(handles.GoodTrials));
      %   MaxFit = Fit_Q(end);
      %   Fit_Q = Fit_Q(P);
    end
    
    
    % ----------------------------------------------------------------------
    % Mark trials as bad based on thresholds.
    
    Trials = (1:nT)';
    %     TrialsLoop = [Trials; Trials(end:-1:1)];
    
    % Trials that were already marked bad and not included.
    BadTrials = setdiff(Trials, handles.GoodTrials);
    nTB = length(BadTrials);
    if nTB ~= nT - nTG
      error('Lost trials.');
    end
    
    TrialCounts.Total(d) = nT;
    TrialCounts.Rejected(d) = nTB;
    
    % Prepare new trial classification but keep old one in case thresholds
    % are changed before saving.
    handles.NewClass = handles.Class;
    % Distance from reference. 'BAD_HeadMotion_Distance'
    % Trial proportion parameter is already taken into account in quantile.
    RejectedTrials = Trials(T_Distance_Q > ...
      handles.DistThresh );
    if ~isempty(RejectedTrials)
      handles.NewClass = AddClass(handles.NewClass, 'BAD_HeadMotion_Distance', RejectedTrials);
    end
    TrialCounts.BadDistance(d) = numel(RejectedTrials);
    
    % Range of motion within trial. 'BAD_HeadMotion_Range'
    RejectedTrials = Trials(T_Range_Q > ...
      handles.DistThresh );
    if ~isempty(RejectedTrials)
      handles.NewClass = AddClass(handles.NewClass, 'BAD_HeadMotion_Range', RejectedTrials);
    end
    TrialCounts.BadMotion(d) = numel(RejectedTrials);
    
    % Fit error. 'BAD_HeadMotion_Fit'
    % Trial proportion parameter is already taken into account in quantile.
    RejectedTrials = Trials(T_Fit_Q > ...
      handles.FitThresh );
    if ~isempty(RejectedTrials)
      handles.NewClass = AddClass(handles.NewClass, 'BAD_HeadMotion_Fit', RejectedTrials);
    end
    TrialCounts.BadFit(d) = numel(RejectedTrials);
    
    
    
    if Verbose
      fprintf('%d: %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.0f%%, %d, %d, %d, %d, %d\n', ...
        d, DistanceQ.Mean(d), DistanceQ.Max(d), DistanceQ.Quantile(d), ...
        MovementQ.Mean(d), MovementQ.Max(d), MovementQ.Quantile(d), MovementQ.Total(d), ...
        FitError.Max(d), TrialCounts.BadDistance(d), TrialCounts.BadMotion(d), ...
        TrialCounts.BadFit(d), TrialCounts.Rejected(d), TrialCounts.Total(d));
    end
    
  end % Dataset loop
  
  if Verbose
    fprintf('\n');
    fprintf('Average (trial counts first converted to proportions):\n');
    fprintf( '  %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.0f%%, %1.2f, %1.2f, %1.2f, %1.2f\n', ...
      mean(DistanceQ.Mean), mean(DistanceQ.Max), mean(DistanceQ.Quantile), ...
      mean(MovementQ.Mean), mean(MovementQ.Max), mean(MovementQ.Quantile), mean(MovementQ.Total), ...
      mean(FitError.Max), mean(TrialCounts.BadDistance ./ TrialCounts.Total), ...
      mean(TrialCounts.BadMotion(d) ./ TrialCounts.Total), ...
      mean(TrialCounts.BadFit ./ TrialCounts.Total), ...
      mean(TrialCounts.Rejected ./ TrialCounts.Total) );
  end
end


