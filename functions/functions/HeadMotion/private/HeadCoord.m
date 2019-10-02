% Head coordinates from distances between pairs of head coils.
function Fids = HeadCoord(CoilDist)
  Fids = zeros(1, 9); % [Na_x, y, z, LE_x, y, z, RE_x, y, z]
  Fids(1) = sqrt( (CoilDist(1)^2 + CoilDist(3)^2 - CoilDist(2)^2/2)/2 );
  Fids(4) = (CoilDist(3)^2 - CoilDist(1)^2) / (4 * Fids(1));
  Fids(5) = sqrt( CoilDist(2)^2/4 - Fids(4)^2 );
  Fids(7:8) = -Fids(4:5);
end
