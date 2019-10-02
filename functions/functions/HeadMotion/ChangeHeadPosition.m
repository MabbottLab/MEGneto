% Matlab version of the CTF command line program changeHeadPos.  Changes
% the dataset .hc and .res4 files such that the initial head coil positions
% is the one provided.  This does not change the data, only the sensor
% coordinates and orientations and the coil coordinates.  Backups of the
% modified files are created (with _Original appended before the file
% extension).

% Fiducials: [Na_x, y, z, LE_x, y, z, RE_x, y, z] coordinates of the head
%            coils in dewar coordinates.

function ChangeHeadPosition(Dataset, Fiducials)
  % Backup files.
  
  % Check that Fiducials are inside helmet, i.e. they are given in dewar
  % coordinates.
  
  % Use dewar sensor positions and orientations.
  
  % Calculate new head coordinates of sensor positions and orientations.
  [Points, Fiducials] = ChangeCoordinates(Points, Fiducials, ...
    Scaling, Orientation, System, Inverse, Vector, Verbose)
  
  % Save new .res4 file.
  % writeRes4(res4File,res4,MAX_COILS)

  % Save new .hc file.
  
  % Any other file?  .infods, etc?
  
end