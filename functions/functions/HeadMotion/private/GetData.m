% Get location data from dataset, in dewar coordinates, in millimeters.
% Get trial classification information from dataset.
function handles = GetData(handles)
  % Set mouse pointer as busy.
  if isfield(handles, 'figure_HeadMotion') && ...
      ishandle(handles.figure_HeadMotion) && handles.figure_HeadMotion ~= 0
    set(handles.figure_HeadMotion, 'Pointer', 'watch');
    drawnow;
    MousePointerChanged = true;
  else
    MousePointerChanged = false;
  end

  % Could use "old way" with Fieldtrip data too, but try to use Fieldtrip
  % user functions instead.
  if handles.Fieldtrip
    ds = ft_read_header(handles.Dataset);
    
    handles.SampleRate = ds.Fs;
    % Number of samples per trial.  Must be constant (as in CTF ds).
    nS = 1 + handles.ft.cfg.trl(1, 2) - handles.ft.cfg.trl(1, 1);
    % ds.nSamples; % would be total number of samples.
    % This is the number of trials we'll get using Fieldtrip to load the data
    % with the desired trial configuration, NOT the trial number in the
    % dataset.
    nT = size(handles.ft.cfg.trl, 1);
    % nCh = ds.nChans;
    % Sensor structure to be corrected later.
    handles.ft.grad = ds.grad;
    % ChannelNames = cellstr(ds.res4.chanNames); % channel names
    if ~isfield(ds.orig.hc, 'dewar')
      % Seems to be a problem with readCTFds in Matlab 2012...
      Fids = openhc(handles.Dataset, true, false)' * 10; % Dewar, not Standard.
      handles.Fiducials(1, :) = Fids(:);
    else
      handles.Fiducials(1, :) = ds.orig.hc.dewar(:) * 10; % [Na_x, y, z, LE_x, y, z, RE_x, y, z] converted from cm to mm.
      % Dewar coords: x points 45 degrees right, y 45 degrees left.
    end
    
    HLC = ds.label(strcmp(ds.chantype, 'headloc'));
    FitChannels = ds.label(strcmp(ds.chantype, 'headloc_gof'));
    %     HLC = {'HLC0011', 'HLC0012', 'HLC0013', 'HLC0021', ...
    %       'HLC0022', 'HLC0023', 'HLC0031', 'HLC0032', 'HLC0033'};
    %     FitChannels = {'HLC0018', 'HLC0028', 'HLC0038'};
    
  else
    
    % Import head localization channels of parsed datasets using CTF m-files.
    ds=readCTFds(handles.Dataset);
    
    handles.SampleRate = ds.res4.sample_rate;
    nS = ds.res4.no_samples;
    nT = ds.res4.no_trials;
    % nCh = ds.res4.no_channels;
    % ChannelNames = cellstr(ds.res4.chanNames); % channel names
    if ~isfield(ds.hc, 'dewar')
      % Seems to be a problem with readCTFds in Matlab 2012...
      Fids = openhc(handles.Dataset, true, false)' * 10;
      handles.Fiducials(1, :) = Fids(:);
    else
      handles.Fiducials(1, :) = ds.hc.dewar(:) * 10; % [Na_x, y, z, LE_x, y, z, RE_x, y, z] converted from cm to mm.
      % Dewar coords: x points 45 degrees right, y 45 degrees left.
    end
    
    HLC = find([ds.res4.senres.sensorTypeIndex] == 13); % ePositionRef = 13
    % eAngleRef = 26; "Orientation of head localization field"
    % eExtractionRef = 27; "Extracted signal from each sensor of field generated by each localization coil"
    FitChannels = find([ds.res4.senres.sensorTypeIndex] == 28); % eFitErr = 28
  end
  
  nC = length(HLC);
  if length(FitChannels) ~= nC/3
    error('Number of fit error channels (%d) doesn''t match number of head coil channels (%d).', length(FitChannels), nC/3);
  end
  if nC == 0
    warning('Continuous head localization may have been off during dataset acquisition, no position channels found.')
  elseif nC < 9
    warning('Found %d position channels, should be 9.', nC);
  elseif nC > 9
    fprintf('Found %d position channels, assuming first 9 are (x,y,z) of the Na, LE and RE coils.', nC);
    HLC(10:end) = [];
    FitChannels(4:end) = [];
    nC = 9;
  end
  
  if handles.Fieldtrip
    cfg = handles.ft.cfg;
    cfg.channel = HLC;
    data = ft_preprocessing(cfg);
    handles.Data = permute(reshape([data.trial{:}] * 1e3, [nC, nS, nT]), [2, 1, 3]); % converted to 'mm'
    % ds.chanunit is 'unknown', but in the current version is in 'm' (2016-03-08).
    cfg.channel = FitChannels;
    data = ft_preprocessing(cfg);
    handles.FitData = permute(reshape([data.trial{:}], [nC/3, nS, nT]), [2, 1, 3]);
    clear cfg data HCL FitChannels
  else
    handles.Data = getCTFdata(ds, [], HLC) * 1e3; % size(Data)=[nS, nC, nT] in dewar coordinates, converted from m to mm.
    handles.FitData = getCTFdata(ds, [], FitChannels); % proportion (no units).
  end
  clear ds
  
  if ndims(handles.Data) < 3 
    if nT ~= 1 || ~all(size(handles.Data) == [nS, nC])
      error('Data dimensions do not match number of samples, channels, trials.');
    end
  elseif ndims(handles.Data) > 3 || ~all(size(handles.Data) == [nS, nC, nT])
    error('Data dimensions do not match number of samples, channels, trials.');
  end
  
  % If head really still this can fail.  Use all trials to find real
  % period.
  handles.HeadSamplePeriod = nS; % Initialized to at least 1 point per trial.
  for t = 1:nT
    % Downsample to real localization sampling rate.  To find it, look for
    % changes in data, but it seems it's common to get repeated values (code
    % must have some condition where it just keeps the same value), so we
    % need to verify carefully.  First, find times of changes.  This already
    % ignores the first few samples until the first change.
    TrueSamples = find(any(diff(handles.Data(:, :, t), 1, 1), 2)) + 1;
    % Then get the time intervals between these changes, and find the
    % smallest "step" between these intervals.  E.g. if we got intervals of
    % 80, 120, 240 samples, the smallest difference between these would be
    % 40, which is the sampling period we were looking for.
    if numel(TrueSamples) <= 1
      continue; % to avoid empty which propagates in min and "erases" our previous good min.
    end
    handles.HeadSamplePeriod = min(handles.HeadSamplePeriod, ...
      min(diff(TrueSamples)));
    %   if numel(handles.HeadSamplePeriod) == 0 % No longer possible.
    %       % Localization sampling period is greater than number of samples per trial.
    %       handles.HeadSamplePeriod = nS;
    %   %else % Do nothing, we have our sampling period.
    %   end
  end
  % Downsample.
  handles.Data = handles.Data(1:handles.HeadSamplePeriod:nS, :, :);
  handles.FitData = handles.FitData(1:handles.HeadSamplePeriod:nS, :, :);
  nS = ceil(nS / handles.HeadSamplePeriod);
  
  % Ignore trials previously classified as bad.
  if handles.Fieldtrip
    % Ignore CTF classification of trials when using Fieldtrip.
    handles.Class = struct('Name', {}, 'Id', {}, 'Count', {}, 'Trials', {});
  else
    handles.Class = opencls(handles.Dataset, false, false); % Not original and quiet.
  end
  handles.GoodTrials = GetGoodTrials(handles.Class, nT);
  

  % Actually, don't do this because it would allow rejecting trials that
  % don't actually exist in the dataset.
  % Create trials for continuous datasets.
  %   if nT == 1
  %     nSNew = floor(1 / handles.SampleRate);
  %     nT = floor(nS / nSNew);
  %     handles.Data = reshape();
  %   end
  
  % Avoid zeros for incomplete trials (acquisition manually interrupted).
  if nT > 1
    nTInc = 0;
    while ~isempty(find(~sum( ...
        handles.Data(:, :, handles.GoodTrials(end - nTInc)), 2 ), 1, 'last')) % Last sample is 0 for all coils.
      nTInc = nTInc + 1; % Reject last trial(s).
    end
    %     Data(:, :, nT+1:end) = [];
    % Mark incomplete trials as bad.
    handles.Class = AddClass(handles.Class, 'BAD_Incomplete', ...
      handles.GoodTrials((end+1-nTInc:end)'));
    handles.GoodTrials(end+1-nTInc:end) = [];
  else % continuous
    % Reject 'zero' samples.  Find last non-zero sample for all coils.
    nS = find(sum( handles.Data, 2 ), 1, 'last'); % Used to be on GoodTrials only, but fails if the single trial has been marked bad.
    handles.Data(nS+1:end, :, :) = [];
    handles.FitData(nS+1:end, :, :) = [];
  end
  
  handles.nS = nS;
  handles.nT = nT;
  
  % Reset mouse pointer.
  pause(0.05); % required so that 'BusyAction' 'cancel' works!
  if MousePointerChanged
    set(handles.figure_HeadMotion, 'Pointer', 'arrow');
  end
  
end