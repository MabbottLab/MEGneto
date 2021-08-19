% Returns the inverse coordinate transformation (Rigid is in head
% coordinates and transformed to non-head coordinates of Locations).
function NewLocations = CorrectRigid(Locations, Rigid)
  switch ndims(Locations)
    case 2 % Single trial dataset.
      NeedReshape = true;
      nS = size(Locations, 1);
      nT = 1;
    case 3
      if size(Locations, 2) == 9
        NeedReshape = true;
        [nS, unused, nT] = size(Locations); % ~
      else % Single trial dataset.
        NeedReshape = false;
        nS = size(Locations, 1);
        nT = 1;
      end
    case 4
      NeedReshape = false;
      [nS, unused, unused2, nT] = size(Locations); % ~, ~
    otherwise
      error('Unexpected dimensions.');
  end
  
  % Reshape.
  Y = reshape(Rigid, [3, 3]); % Points (coils) are columns.

  NewLocations = zeros(nS, 9, nT);
  for t = 1:nT
    for s = 1:nS
      if NeedReshape
        [XO, XR] = RigidCoordinates(reshape(Locations(s, :, t), [3, 3]));
      else
        [XO, XR] = RigidCoordinates( squeeze(Locations(s, :, :, t)) );
      end
      % Inverse transformation going from Y head coordinates to X dewar.
      NewLocations(s, :, t) = reshape(XR * Y + XO(:, [1, 1, 1]), [1, 9, 1]);
    end
  end
end
