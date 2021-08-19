% Max of the 3 coil norms.
function MaxNorms = MaxCoilNorms(Locations)
  Norms = [sqrt(sum(Locations(:, 1:3, :).^2, 2)), ...
    sqrt(sum(Locations(:, 4:6, :).^2, 2)), ...
    sqrt(sum(Locations(:, 7:9, :).^2, 2))];
  
  MaxNorms = max(Norms, [], 2);
end