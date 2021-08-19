% Overall geometric median location of each coil.
function MedianLoc = MedianLocation(Data, GoodTrials)
  nSxGT = size(Data, 1) * length(GoodTrials);
  if nSxGT == 0 % If all trials rejected.
    MedianLoc = zeros(1, 9);
  else
    MedianLoc = reshape( GeoMedian( ...
      reshape(permute( Data(:, :, GoodTrials), ...
      [1, 3, 2]), [nSxGT, 3, 3]), 1e-3 ) , [1, 9]);
  end
end
