function [Points, Fiducials] = ChangeCoordinates(Points, Fiducials, ...
    Scaling, Orientation, System, Inverse, Vector, Verbose)
  % Change to/from CTF or Curry coordinates given MEG fiducial coordinates.
  %
  % [Points, Fiducials] = ChangeCoordinates(Points, Fiducials, ...
  %     Scaling, Orientation, System, Inverse, Vector, Verbose)
  % [Origin, Rotation] = ChangeCoordinates([], Fiducials, ...)
  %
  % Outputs:
  % Normally returns the Points and Fiducials in the new coordinates
  % system. However, if Points input is empty, returns [Origin, Rotation]
  % instead, in scaled coordinates.
  %
  % Inputs:
  % Only Fiducials is required, all other arguments are optional.
  % Points = Array of 3d points (or vectors see Vector option below) to
  %   be transformed (nV x 3).
  % Fiducials = [Na_x, Na_y, Na_z; LE_x, LE_y, LE_z; RE_x, RE_y, RE_z];
  %   Points and Fiducials must have same units since we translate.
  % Scaling can be a scalar or a 1x3 vector (e.g. voxel size).  It is the
  %   factor(s) by which to multiply to get cm.  So e.g. if input Points
  %   are in mm, Scaling should be 0.1.
  % Orientation is +1 or -1 indicating if the new and old coordinate
  %   systems have the same or opposite handedness respectively. Since both
  %   possible new systems are right-handed (see System below), -1
  %   indicates that the old system was left-handed.
  % System should be 'CTF' (ALS) (default) or 'Curry' (LPS) coordinates.
  % Inverse is a boolean that indicates if we want the forward (default) or
  %   inverse coordinate transformation applied to the Points.
  % Vector is a boolean that indicates if the points in Points are actually
  %   vectors (e.g. orientation vectors) and should not be translated, only
  %   scaled and rotated.
  %
  % More about Orientation:
  % Orientation here is used in the mathematical sense: the handedness of a
  % coordinate system. It is necessary to indicate when it changes because
  % the fiducials only define a plane and there is no way to know which
  % perpendicular direction is superior: any coordinate system is simply a
  % vector space and appears to have positive orientation (i x j = +k) on
  % its own. So if we know that the original coordinate system was
  % left-handed, we know that Z should be -X x Y.
  %
  % Marc Lalancette, The Hospital for Sick Children, Toronto, Canada.
  % 2014-02-06

  % Parse inputs and assign default values.
  if ~exist('Verbose', 'var') || isempty(Verbose)
    Verbose = 1;
  end
  if ~exist('Fiducials', 'var') || isempty(Fiducials)
    error('Fiducials (second input) required.');
  end
  if ~ismember(size(Points, 2), [0, 3])
    error('Points should be an N by 3 matrix (or empty).');
  end
  % Removed as unnecessary input:
  %   if ~exist('nV', 'var') || isempty(nV)
  nV = size(Points, 1);
  %   end
  if ~exist('Orientation', 'var') || isempty(Orientation) ...
      || Orientation == 0 % If misinterpreted as a boolean.
    Orientation = 1;
  elseif ~ismember(Orientation, [1, -1])
    error('Orientation should be +1 or -1 indicating a right or left-handed original coordinate system.');
  end
  if ~exist('Scaling', 'var') || isempty(Scaling)
    Scaling = ones(1, 3);
  elseif numel(Scaling) == 1
    % Apply scaling uniformly.
    Scaling = Scaling(ones(1, 3));
  elseif ~all(size(Scaling) == [1, 3])
    error('Scaling must be a scalar or a 1 by 3 vector. size(Scaling) = %s', ...
      num2str(size(Scaling)) )
  end
  if ~exist('System', 'var') || isempty(System)
    System = 'CTF';
  end
  if ~exist('Inverse', 'var') || isempty(Inverse)
    Inverse = false;
  end
  if ~exist('Vector', 'var') || isempty(Vector)
    Vector = false;
  end
  
  
  % Resize first since voxels can be anisotropic (non-cubic).
  % Have to apply the same transformation to fiducials and Points (below).
  Fiducials = Fiducials .* Scaling(ones(1, 3), :);
  
  % Prepare coordinate system.
  if strcmpi(System, 'CTF')
    % Source: getAffineVox2CTF.m  Had bug where Y and Z were not unit length.
    Origin = (Fiducials(2, :) + Fiducials(3, :))/2;
    X = Fiducials(1, :) - Origin;
    X = X/norm(X);
    Y = Fiducials(2, :) - Origin; % Not yet perpendicular to X in general.
    Y = Y/norm(Y);
    Z = cross(X, Y);
    Z = Z/norm(Z); % Necessary.
    Y = cross(Z, X); % Doesn't go through ears anymore in general.
    % Should be ok, but check.
    if abs(norm(Y) - 1) > 1e-15
      error('X and Z aren''t perpendicular, norm(Y) - 1 = %g.', norm(Y) - 1)
    end
  elseif strcmpi(System, 'Curry')
    X = Fiducials(2, :) - Fiducials(3, :);
    X = X/norm(X);
    % Origin is on ears line, where a perpendicular drops from the nasion.
    Origin = Fiducials(2, :) + (Fiducials(1, :) - Fiducials(2, :)) * X' * X;
    Y = -Fiducials(1, :) + Origin;
    Y = Y/norm(Y);
    Z = cross(X, Y);
    % Should be ok, but check.
    if abs(norm(Z) - 1) > 1e-15
      error('X and Y aren''t perpendicular, norm(Z) - 1 = %g.', norm(Z) - 1)
    end
  else
    error('Unrecognized coordinate system: %s.', System)
  end
  
  % If the original system orientation was negative (left-handed),
  % Z = X x Y would also produce a left-handed system and the Z axis would
  % point in the inferior direction. So we must do:
  Z = Orientation * Z;
  Rotation = [X; Y; Z];
  
  if Vector
    % Points are actually vectors (e.g. orientation vectors) so do not
    % translate.
    Origin = zeros(1, 3);
  end
  
  if isempty(Points)
    % Return [Origin, Rotation] instead, but in scaled
    % coordinates so scaling must be performed also.
    Points = Origin;
    Fiducials = Rotation;
    if ~all(Scaling == ones(1, 3)) && Verbose
      fprintf('Scaling (%g %g %g) MUST be performed on Points before applying translation and rotation.\n', ...
        Scaling)
    end
  else
    % Return fiducials in new coordinate system (regardless of Inverse).
    % Fiducials already scaled.
    Fiducials = (Fiducials - Origin(ones(1, 3), :)) * Rotation';
    if ~Inverse
      % Apply forward coordinate transformation: scale, translate, rotate.
      Points = ( Points .* Scaling(ones(1, nV), :) ...
        - Origin(ones(1, nV), :) ) * Rotation';
    else
      % Apply inverse coordinate transformation: rotate', -translate, /scale.
      Points = ( Points * Rotation + Origin(ones(1, nV), :) ) ...
        ./ Scaling(ones(1, nV), :);
    end
  end
  
end