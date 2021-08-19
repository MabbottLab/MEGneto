% Copy of part of the external ChangeCoordinates function, simplified and
% optimized for speed given it gets called many times.  Returns the origin
% and rotation matrix.
function [O, R] = RigidCoordinates(FidsColumns)
  R = zeros(3);
  O = (FidsColumns(:, 2) + FidsColumns(:, 3))/2;
  R(:, 1:2) = FidsColumns(:, 1:2) - O(:, [1, 1]);
  %R(:, 3) = cross(R(:, 1), R(:, 2));
  R(:, 3) = [R(2, 1)*R(3, 2) - R(3, 1)*R(2, 2), -R(1, 1)*R(3, 2) + R(3, 1)*R(1, 2), R(1, 1)*R(2, 2) - R(2, 1)*R(1, 2)];
  %R(:, 2) = cross(R(:, 3), R(:, 1));
  R(:, 2) = [R(2, 3)*R(3, 1) - R(3, 3)*R(2, 1), -R(1, 3)*R(3, 1) + R(3, 3)*R(1, 1), R(1, 3)*R(2, 1) - R(2, 3)*R(1, 1)];
  % Normalize x, y, z.
  R = bsxfun(@rdivide, R, sqrt(sum(R.^2, 1)));
end
