% Maximum distance within a spherical volume: translation + rotation about
% origin.  The formula used is equivalent to 2*R*sin(a/2) where a is the
% angle of the single rotation equivalent to the original
% translation+rotation, and R is the distance from that axis to the
% furthest edge of the sphere.  For the sphere radius, we use the maximum
% coil distance from the origin in the rigid body reference.  This is done
% independently for each time sample.
function D = RigidDistances(Locations, Reference)
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
  Y = reshape(Reference, [3, 3]); % Points (coils) are columns.
  
  % Reference "head origin" and inverse "orientation matrix".
  [YO, YR] = RigidCoordinates(Y);
  % Sphere radius.
  r = max( sqrt(sum((Y - YO(:, [1, 1, 1])).^2, 1)) ); 
  if any(YR(:)) % any ignores NaN and returns false for empty.
    YI = inv(YR); % Faster to calculate inverse once here than "/" in loop.
  else
    YI = YR;
  end
  
  %   SinHalf = zeros([nS, 1, nT]);
  %   Axis = zeros([nS, 3, nT]);
  D = zeros([nS, 1, nT]);
  for t = 1:nT
    for s = 1:nS
      if NeedReshape
        [XO, XR] = RigidCoordinates(reshape(Locations(s, :, t), [3, 3]));
      else
        [XO, XR] = RigidCoordinates( squeeze(Locations(s, :, :, t)) );
      end
      % Translation from X "head origin" to Y "head origin".
      T = XO - YO;
      % Rotation from X to Y (both with their "head origin" subtracted, so
      % it is a rotation around an axis through the real origin).
      R = XR * YI; % %#ok<MINV>
      %       % Correct for numerical errors.  Not necessary now that we
      %       % use w for SinHalf.
      %       R(:, 3) = cross(R(:, 1), R(:, 2));
      %       R(:, 2) = cross(R(:, 3), R(:, 1));
      %       R = bsxfun(@rdivide, R, sqrt(sum(R.^2, 1)));
      
      TrR = trace(R);
      %       % Sine of half the rotation angle.
      %       SinHalf = sqrt(3 - TrR) / 2;
      %       % For very small angles, this formula is not accurate compared to w,
      %       % since diagonal elements are around 1, and eps(1) = 2.2e-16.  This
      %       % will be the order of magnitude of non-diag. elements due to errors.
      %       % So we should get SinHalf from w.
      % Rotation axis with amplitude = SinHalf (like in rotation quaternions).
      w = [R(3, 2) - R(2, 3); R(1, 3) - R(3, 1); R(2, 1) - R(1, 2)] / ...
        (2 * sqrt(1 + TrR));
      SinHalf = sqrt(sum(w.^2));
      TNormSq = sum(T.^2);
      % Maximum sphere distance for translation + rotation, as described
      % above.
      D(s, 1, t) = sqrt( TNormSq + (2 * r * SinHalf)^2 + ...
        4 * r * sqrt(TNormSq * SinHalf^2 - (T' * w)^2) );
      % CHECK should be comparable AND >= to max coil movement.
    end
  end
    
end