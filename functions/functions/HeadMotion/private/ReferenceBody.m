% Calculate a reference "rigid body" for the relative positions of the head
% coils, to deal with head coils possibly falling off, or moving, and to
% later calculate movement amplitude based on whole head volume, not simply
% at coil locations.  Also plot inter-coil distance deviation from start
% and stable interval boundaries, and display some related values.
function handles = ReferenceBody(handles)
  % Set mouse pointer as busy.
  if isfield(handles, 'figure_HeadMotion') && ...
      ishandle(handles.figure_HeadMotion) && handles.figure_HeadMotion ~= 0
    set(handles.figure_HeadMotion, 'Pointer', 'watch');
    drawnow;
    GUI = true;
  else
    GUI = false;
  end
  
  % Changed to use all trials instead of "good trials", because this is for
  % finding a good model of the coil positions before they move on the
  % skin.  Even if we reject the first trials, we should trust the coil
  % positions relative to each other more at the start.
  nSxT = handles.nS * handles.nT;
  %   nSxGT = handles.nS * length(handles.GoodTrials);
  % Calculate distances between pairs of head coils.
  CoilDiffs = [handles.Data(:, 1:3, :) - handles.Data(:, 4:6, :), ...
    handles.Data(:, 4:6, :) - handles.Data(:, 7:9, :), ...
    handles.Data(:, 7:9, :) - handles.Data(:, 1:3, :)];
  CoilDistances = [sqrt(sum(CoilDiffs(:, 1:3, :).^2, 2)), ...
    sqrt(sum(CoilDiffs(:, 4:6, :).^2, 2)), ...
    sqrt(sum(CoilDiffs(:, 7:9, :).^2, 2))];
  CoilDiffs = reshape(permute(CoilDiffs, [1, 3, 2]), [nSxT, 9]); % [nSxT, 9]
  % For plotting.
  %   CoilDiffs = bsxfun(@minus, CoilDiffs, CoilDiffs(1, :));
  
  %   CoilDistances = CoilDistances(:, :, handles.GoodTrials);
  CoilDistances = reshape(permute(CoilDistances, [1, 3, 2]), [nSxT, 3]); % [nSxT, 3]
  % Distance deviations compared to initial distances.  Used for plotting
  % and for simple method of finding initial distances.
  InterDDev = bsxfun(@minus, CoilDistances, CoilDistances(1, :));
  
  % Calculate uncertainty.  Use maximum-likelihood estimate of variance,
  % ignoring zeros, since it looks like the distribution is Maxwell (chi) +
  % a large number of zeros.  The algorithm probably does something like
  % "if small difference, keep same position", which generates these
  % duplicate values.
  InterDAdjDiff = diff(CoilDiffs, 1, 1);
  InterDAdjDiff = [sqrt(sum(InterDAdjDiff(:, 1:3, :).^2, 2)), ...
    sqrt(sum(InterDAdjDiff(:, 4:6, :).^2, 2)), ...
    sqrt(sum(InterDAdjDiff(:, 7:9, :).^2, 2))];
  NonZeroCounts = sum(InterDAdjDiff > 0, 1);
  Sigmas = sqrt(sum(InterDAdjDiff.^2, 1) ./ NonZeroCounts); % 1/3 inside the sqrt would give for each component.  This is for the distances.
  
  % Use inter-quartile range of adjacent sample inter-coil distance
  % distributions.  Dosen't see as many changes as "distance-diff", but
  % then it's positive only, so the relation between IQR and STD is not as
  % simple, which is why we went with an estimate of the variance based on
  % maximum-likelihood estimate on the Maxwell (chi) distribution.
  %   InterDAdjDiff2 = diff(InterD, 1, 1);
  %   IQR = InterQuartileRange(InterDAdjDiff);
  
  % Check if max variation is big enough to investigate further and
  % possibly restrict the range to build reference rigid body.  Use 3 sigma
  % based on max inter-coil IQR.
  
  % [NO LONGER: Use a threshold one order of magnitude smaller than the
  % desired movement threshold for the head.]
  %   MoveThresh = str2double(get(handles.edit_Thresh, 'String')) / 10;
  %   MoveThresh = 1.57 * max(IQR); % 3 sigma of parent Gaussian = 1.57 IQR of difference distribution (convolution of 2 Gaussians = Gaussian with sqrt(2)sigma std).
  %   MoveThresh = 3 * max(Sigmas);
  if GUI
    MoveThresh = str2double(get(handles.edit_InterCoilThresh, 'String'));
  else
    MoveThresh = handles.DistThresh / 10;
  end
  
  % Use simple threshold way to select the start of the dataset where the
  % inter-coil distances are "stable".  Otherwise, uses a partitioning
  % algorithm, but it turns out that the coils move much more than
  % anticipated and the partitioning algorithm is not very useful, on top
  % of not being very efficient on older computers.  So it was commented
  % out in the code since it also removes the requirement for mex files.
  %   SimpleRigid = true;
  %   if SimpleRigid
  
  % Keep the same "interval" organization for now, but with just 1 or 2.
  FirstMoveSample = find(any(abs(InterDDev) > MoveThresh, 2), 1, 'first');
  if ~isempty(FirstMoveSample)
    I = [1, FirstMoveSample - 1; FirstMoveSample, nSxT];
    nI = 2;
  else
    I = [1, nSxT];
    nI = 1;
  end
  FirstMove = 2; % Index of interval where the move occured.
  
  %   else % Partitioning method, requires mex files.
  %   if any(max(CoilDistances, [], 1) - min(CoilDistances, [], 1) > MoveThresh)
  %     % Partition into stationary intervals (external function).
  %     % First, remove duplicate adjacent values, which account for a quarter
  %     % of the samples (duplicate across the 3 coils at the same time)!
  %     NonDuplicates = [true; any(diff(CoilDistances, 1, 1) ~= 0, 2)];
  %     NonDupIx = find(NonDuplicates);
  %     % Use a second for the minimum "stable interval" time. (30)
  %     nMin = min(nSxT, round(2 * handles.SampleRate / handles.HeadSamplePeriod));
  %     % Removing duplicate values didn't help much, could also "dither" for
  %     % repeated values in one coil, but not much improvement either and
  %     % makes results non-deterministic.  So instead, use a smaller p-value
  %     % and larger nMin.
  %     [I, nI] = PartitionTimeSeries(CoilDistances(NonDupIx, :), nMin, 0.0001); % 6 sec for 18k points, 10 minutes of data.
  %     %  + min(Sigmas)/2 * rand(length(NonDupIx), 3)
  %     I(:, 1) = NonDupIx(I(:, 1));
  %     I(:, 2) = [I(2:end, 1) - 1; nSxT];
  %
  %     % Find first coil move above threshold.
  %     InterM = zeros(nI, 3);
  %     for i = 1:nI
  %       InterM(i, :) = median(CoilDistances(I(i, 1):I(i, 2), :), 1);
  %     end
  %     % Index of interval where the move occured (at the start of interval).
  %     FirstMove = find(any( ...
  %       abs(bsxfun(@minus, InterM, InterM(1, :))) > MoveThresh, ...
  %       2), 1, 'first'); % index into nI
  %     % else FirstMove variable does not exist, this is used as check below.
  %   else
  %     I = [1, nSxT];
  %     nI = 1;
  %   end
  %   if ~exist('FirstMove', 'var') || isempty(FirstMove)
  %     FirstMove = nI + 1;
  %   end
  %   end % Simple or Partition
  
  % Build the rigid body.  Use median inter-coil distances up to the first
  % move above threshold.
  RigidD = median(CoilDistances(I(1, 1):I(FirstMove-1, 2), :), 1);
  % Get rigid body coil coordinates (in head coordinates now, dewar later).
  handles.Rigid = HeadCoord(RigidD);
  
  
  % ----------------------------------------------------------------------
  % Fill in inter-coil distance plot and related text boxes.
  if GUI
    
    Darker = 153/255;
    Dark = 204/255;
    Blue = [0, 0, Dark];
    Red = [Dark, 0, 0];
    Purple = [Darker, 0, Dark];
    Green = [0, Darker, 0];
    Orange = [Dark, Dark/2, 0];
    Colors = {Blue, Green, Purple};
    Padding = 1.1;
    
    X = 1 + (0:nSxT-1)'/handles.nS;
    
    % Get trial and sample, and x axis coordinate for interval starts.
    IT = floor(I(:, 1) / handles.nS) + 1; %handles.GoodTrials()
    IS = mod(I(:, 1), handles.nS);
    IX = (IT) + (IS - 1) / handles.nS;
    if FirstMove > nI
      FirstMoveX = -1;
    else
      FirstMoveX = IX(FirstMove);
    end
    IX = reshape(IX(:, [1, 1])', [], 1); % Duplicate for zig-zag plot.
    % y coordinates are below and above y axis limits, in a zig-zag way:
    % down, up, up, down, etc.
    %   set(handles.axes_InterCoil, 'YLimMode', 'Manual');
    %   YLim = get(handles.axes_InterCoil, 'YLim') + [-1, 1];
    YLim = Padding * [min(InterDDev(:)), max(InterDDev(:))];
    set(handles.axes_InterCoil, 'YLim', YLim);
    IY = Padding * [YLim'; YLim(2); YLim(1)];
    IY = IY(:, ones(ceil(nI/2), 1));
    IY = IY(1:(2*nI));
    
    if ~isfield(handles, 'XLabelInterCoil')
      set(gcf, 'CurrentAxes', handles.axes_InterCoil);
      AxPos = get(handles.axes_InterCoil, 'Position');
      handles.XLabelInterCoil = text(AxPos(3)/2, 10 - handles.TextPos, ...
        'Trials', handles.TextOptions{:}, ...
        'Parent', handles.axes_InterCoil, 'VerticalAlignment', 'bottom'); %, 'Visible', 'off'
      handles.YLabelInterCoil = text(-handles.TextPos, AxPos(4)/2, ...
        'Inter-coil distance deviations (mm)', handles.TextOptions{:}, ...
        'Parent', handles.axes_InterCoil, 'Rotation', 90, ...
        'VerticalAlignment', 'top'); %, 'Visible', 'off'
      
      handles.Plot_InterCoil = zeros(3, 1);
      for Coil = 1:3
        %       for Coord = 1:3
        %         handles.Plot_InterCoil((Coil-1)*3 + Coord) = ...
        %           plot(X, CoilDiffs(:, (Coil-1)*3 + Coord), 'Color', Colors{Coil});
        %       end
        handles.Plot_InterCoil(Coil) = ...
          plot(X, InterDDev(:, Coil), 'Color', Colors{Coil});
      end
      
      handles.Plot_InterCoil_Zero = plot([X(1), X(end)], [0, 0], ':k');
      
      % Draw vertical lines for boundaries of stable intervals.
      handles.Plot_I = plot(IX, IY, '-', 'Color', Orange);
      handles.Plot_FirstMove = plot(FirstMoveX * [1, 1], Padding * YLim, ...
        '--', 'LineWidth', 2, 'Color', Red);
    else % Plot exists, update data only.
      for Coil = 1:3
        set(handles.Plot_InterCoil(Coil), 'XData', X, ...
          'YData', InterDDev(:, Coil));
      end
      set(handles.Plot_InterCoil_Zero, 'XData', [X(1), X(end)]);
      set(handles.Plot_I, 'XData', IX, 'YData', IY);
      set(handles.Plot_FirstMove, 'XData', FirstMoveX * [1, 1], ...
        'YData', Padding * YLim);
    end
    
    % Hide plot if it's not already showing (typically when loading a
    % dataset).
    if ~get(handles.togglebutton_InterCoil, 'Value')
      set(handles.axes_InterCoil, 'Visible', 'off');
      set(get(handles.axes_InterCoil, 'children'), 'Visible', 'off');
    end
    % Set proper X limits based on slider values and make background
    % transparent to see trial highlighting.
    XLim = get(handles.slider_Trials, 'Value') + [0, ...
      handles.nT - get(handles.slider_Trials, 'Max') + 1];
    set(handles.axes_InterCoil, 'XLim', XLim, 'Color', 'none');
    
    % Fill text boxes.
    
    set(handles.edit_NALE, 'String', num2str(Sigmas(1), handles.NumFormatSmall));
    set(handles.edit_LERE, 'String', num2str(Sigmas(2), handles.NumFormatSmall));
    set(handles.edit_RENA, 'String', num2str(Sigmas(3), handles.NumFormatSmall));
    % Factor of 0.5242 converts IQR to std of non-diff distribution, assuming
    % Gaussian.
    %   set(handles.edit_NALE, 'String', num2str(0.5242 * IQR(1), handles.NumFormatSmall));
    %   set(handles.edit_LERE, 'String', num2str(0.5242 * IQR(2), handles.NumFormatSmall));
    %   set(handles.edit_RENA, 'String', num2str(0.5242 * IQR(3), handles.NumFormatSmall));
    
    set(handles.edit_InterCoilThresh, 'String', ...
      num2str(MoveThresh, handles.NumFormat));
    set(handles.edit_nI, 'String', ...
      num2str(nI, handles.NumFormatInt));
    if FirstMove <= nI
      set(handles.togglebutton_InterCoil, 'ForegroundColor', Red, ...
        'FontWeight', 'Bold');
      set(handles.edit_InterCoilTrial, 'String', ...
        num2str(IT(FirstMove), handles.NumFormatInt));
      set(handles.edit_InterCoilSample, 'String', ...
        [num2str(IS(FirstMove), handles.NumFormatInt), ' of ', ...
        num2str(handles.nS, handles.NumFormatInt)]);
    else
      if nI > 1
        set(handles.togglebutton_InterCoil, 'ForegroundColor', Orange, ...
          'FontWeight', 'Normal');
      else
        set(handles.togglebutton_InterCoil, 'ForegroundColor', [0, 0, 0], ...
          'FontWeight', 'Normal');
      end
      set(handles.edit_InterCoilTrial, 'String', '-');
      set(handles.edit_InterCoilSample, 'String', '-');
    end
    
    % Reset mouse pointer.
    pause(0.05); % required so that 'BusyAction' 'cancel' works!
    set(handles.figure_HeadMotion, 'Pointer', 'arrow');
  end
  
end
