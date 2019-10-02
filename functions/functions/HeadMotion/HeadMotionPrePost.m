function [HeadMov, CoilMov] = HeadMotionPrePost(Dataset, Verbose)
  % Returns the head movement in cm, from pre/post localizations in a MEG dataset.
  %
  % [HeadMov, CoilMov] = HeadMotionPrePost(Dataset, Verbose)
  %
  % Head movement is calculated from the "head zeroing" datasets collected
  % before and after the main MEG recording: hz.ds, hz2.ds.  The movement
  % is based on head coil positions, but is calculated to represents the
  % maximal motion of any point on the head, which is modeled as a sphere
  % for this purpose.
  %
  % CoilMov is the maximal difference of head coil pair distances between
  % pre and post recordings.  It serves as a check for the coils having
  % moved themselves.  It is therefore highly recommended to verify this
  % value before interpreting the head movement.
  %
  % Verpose [default true]: Optional argument, if true, displays a little
  % more feedback during execution.
  %
  % Marc Lalancette, The Hospital for Sick Children, Toronto, Canada.
  % 2015-01-09
  
  % .hc file is either pre post or mean, selected by user (in rp file
  % likely).
  
  % .rp file used for "head zeroing" needed, but it is the same as the .acq file
  % in the hz.ds sub-datasets.
  
  if nargin < 2
    Verbose = true;
  end
  
  if Verbose
    fprintf('This command should not take more than a few seconds.  \nIf it appears frozen, you will want to "kill" it with ctrl-c here \n(or the kill command in a terminal).\n');
    drawnow;
  end
  
  % This would work but updates the dataset, so we'd need to do it 3 times,
  % the last time with -overall.
  %     Command = sprintf('calcHeadPos -hzDs %s %s %s', ...
  %       fullfile(Dataset, [hz, '.ds']), Dataset, ...
  %       fullfile(Dataset, [hz, '.ds'], [hz, '.acq']));
  %     [Status, Message] = system(Command, '-echo');
  
  % "Nicer" to run this and parse the output.  We pass 0 to the command
  % which is the option to exit without modifying the dataset.
  Command = sprintf('calcHeadPos %s %s', Dataset, ...
    fullfile(Dataset, 'hz.ds', 'hz.acq'));
  if Verbose > 1
    [Status, Message] = system(['echo "0" | ', Command], '-echo');
  else
    [Status, Message] = system(['echo "0" | ', Command]); % , '-echo'
  end
  if ~Status && Verbose
    fprintf('calcHeadPos command completed succesfully.\n');
  end
  
  % Convert the single line string to cell array.
  MessageC = textscan(Message, '%s', 'delimiter', sprintf('\n'));
  MessageC = MessageC{1};
  
  % Extract "means" coordinates for pre and post.
  Coord = zeros(2, 3, 3); % 2 Hz, 3 dimensions, 3 coils.  Coils are "columns".
  Coils = {'nasion', 'left ear', 'right ear'};
  hz = {'hz.ds', 'hz2.ds'};
  for h = 1:2
    iHz = find(~cellfun(@isempty, strfind(MessageC, hz{h})), 1);
    % hz.ds can have one or multiple trials.  If multiple, the calcHeadPos
    % shows the position for each trial then "Means".  Otherwise, only
    % "Trial 1".
    iMeans = find(~cellfun(@isempty, strfind(MessageC(iHz+1:end), 'Means')), 1);
    if isempty(iMeans)
      iMeans = find(~cellfun(@isempty, strfind(MessageC(iHz+1:end), 'Trial 1')), 1);
    end
    for c = 1:3
      iCoil = find(~cellfun(@isempty, strfind(MessageC(iHz+iMeans+1:end), Coils{c})), 1);
      for i = 1:3
        x = textscan(MessageC{iHz+iMeans+iCoil+i}, '%*s %*s %f');
        if isempty(x)
          error('Didn''t find coordinate %d, coil %d, hz %d.', i, c, h);
        end
        Coord(h, i, c) = x{1};
      end
    end
  end
  
  
  % Dewar coordinates, in cm.  HeadMotionTool uses mm, but shouldn't matter
  % for RigidDistances.
  % Important that the first dim be singleton, would be samples in
  % HeadMotionTool.
  HeadMov = RigidDistances(Coord(2, :, :), Coord(1, :, :));
  
  if nargout > 1
    CoilDiffs = [Coord(:, :, 1) - Coord(:, :, 2), ...
      Coord(:, :, 2) - Coord(:, :, 3), ...
      Coord(:, :, 3) - Coord(:, :, 1)];
    CoilDistances = [sqrt(sum(CoilDiffs(:, 1:3).^2, 2)), ...
      sqrt(sum(CoilDiffs(:, 4:6).^2, 2)), ...
      sqrt(sum(CoilDiffs(:, 7:9).^2, 2))];
    
    CoilDistDev = abs(CoilDistances(2, :) - CoilDistances(1, :));
    CoilMov = max(CoilDistDev);
  end
end


