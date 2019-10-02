function M = GeoMedian(X, Precision)
  
  % M = GeoMedian(X, Precision)
  %
  % Calculate the geometric median: the point that minimizes sum of
  % Euclidean distances to all points.  size(X) = [n, d, ...], where n is
  % the number of data points, d is the number of components for each point
  % and any additional array dimension is treated as independent sets of
  % data and a median is calculated for each element along those dimensions
  % sequentially; size(M) = [1, d, ...].  This is an approximate iterative
  % procedure that stops once the desired precision is achieved.  If
  % Precision is not provided, 1e-4 of the max distance from the centroid
  % is used.
  % 
  % Weiszfeld's algorithm is used, which is a subgradient algorithm; with
  % (Verdi & Zhang 2001)'s modification to avoid non-optimal fixed points
  % (if at any iteration the approximation of M equals a data point).
  %
  % Marc Lalancette, last updated 2012-05

  nDims = ndims(X);
  XSize = size(X);
  n = XSize(1);
  d = XSize(2);
  if nDims > 3
    nSets = prod(XSize(3:nDims));
    X = reshape(X, [n, d, prod(XSize(3:nDims))]);
  elseif nDims == 3
    nSets = XSize(3);
  else
    nSets = 1;
  end
  
  % For better stability, center and normalize the data.
  Centroid = mean(X, 1);
  Scale = max(max(abs(X), [], 1), [], 2); % [1, 1, nSets]
  X = bsxfun(@rdivide, bsxfun(@minus, X, Centroid), Scale); % (X - Centroid(ones(n, 1), :, :)) ./ Scale(ones(n, 1), ones(d, 1), :);
  
  if ~exist('Precision', 'var') || isempty(Precision)
    Precision = 1e-4 * ones(1, 1, nSets);
  else
    Precision = bsxfun(@rdivide, Precision, Scale); % Precision ./ Scale; % [1, 1, nSets]
  end
  
  % Initial estimate: median in each dimension separately.  Though this
  % gives a chance of picking one of the data points, which requires
  % special treatment.
  M2 = median(X, 1);
  
  % It might be better to calculate separately each independent set,
  % otherwise, they are all iterated until the worst case converges.
  for s = 1:nSets
    
    % For convenience, pick another point far enough so the loop will always
    % start.
    M = bsxfun(@plus, M2(:, :, s), Precision(:, :, s));
    % Iterate.
    while  sum((M - M2(:, :, s)).^2 , 2) > Precision(s)^2  % any()scalar
      M = M2(:, :, s); % [n, d]
      % Distances from M.
      %       R = sqrt(sum( (M(ones(n, 1), :) - X(:, :, s)).^2 , 2 )); % [n, 1]
      R = sqrt(sum( bsxfun(@minus, M, X(:, :, s)).^2 , 2 )); % [n, 1]
      % Find data points not equal to M, that we use in the computation
      % below.
      Good = logical(R);
      nG = sum(Good);
      if nG % > 0
        %       D = sum( (M(ones(nG, 1), :) - X(Good, :, s)) ./ R(Good, ones(d, 1)) , 1 ); % [1, d, 1]
        D = sum( bsxfun(@rdivide, bsxfun(@minus, M, X(Good, :, s)), R(Good)) , 1 ); % [1, d, 1]
        %       DNorm = sqrt(sum( D.^2 , 2 )); % scalar
        %       W = sum(1 ./ R, 1); % scalar. Sum of "weights" (in one viewpoint of this problem).
      else % all points are in the same location 
        % Above formula would give error due to second bsxfun on empty.
        D = 0;
      end
      
      % New estimate. 
      % Note the possibility of D = 0 and (n - nG) = 0, in which case 0/0
      % should be 0, but here gives NaN, which the max function ignores,
      % returning 0 instead of 1. This is fine however since this
      % multiplies D (=0 in that case).
      M2(:, :, s) = M - max(0, 1 - (n - nG)/sqrt(sum( D.^2 , 2 ))) * ...
        D / sum(1 ./ R, 1);
    end
    
  end
  
  % Go back to original space and shape.
  %   M = M2 .* Scale(1, ones(d, 1), :) + Centroid;
  M = bsxfun(@times, M2, Scale) + Centroid;
  if nDims > 3
    M = reshape(M, [1, XSize(2:end)]);
  end
  
end
  