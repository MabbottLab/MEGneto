function y = d2b(x, nBits)
  
  % Convert decimal numbers into a logical array, least significant digits
  % first. The function is meant to work on integers and the fractional
  % part is truncated. Adds a trailing dimension for binary digits.
  % Optional argument nBits pads or truncates most significant bits so that
  % the size of the added dimension is nBits.
  %
  % d2b by Zacharias Voulgaris, from Matlab File Exchange.
  % Modified 2012-04 by Marc Lalancette to work on arrays.
  
  %   % except for row vector where bits are along columns. (NO LONGER)
  
  % Number of binary digits.
  if ~exist('nBits', 'var') || isempty(nBits)
    nBits = floor( log2(max(x(:))) ) + 1;
  end
  
  % Find which dimension to use.
  % There is no second dimension if column vector!
  nD = ndims(x); % Everything is at least a 2-dim array to Matlab, and additional trailing singleton dimensions are ignored.
  InputSize = size(x);
  if nD == 2 && InputSize(2) == 1
    nD = nD - 1;
    InputSize(2) = [];
  end
  
  %   % Use first dimension if row vector.
  %   if nD == 2 && InputSize(1) == 1
  %     x = x';
  %     InputSize(1) = [];
  %     nD = nD - 1;
  %     Transposed = true;
  %   else
  %     Transposed = false;
  %   end
  
  % Initialize output array.
  y = false([InputSize, nBits]);
  InputDims(1:nD) = {':'};

  % In case the inputs are not integers.
  x = floor(x);

  % Compute the bits.
  for b = 1:nBits
    x = x ./ 2;
    y(InputDims{:}, b) = x ~= fix(x);
    x = floor(x);
  end
  % This was relatively equivalent in terms of speed.
  %     for b = 1:nBits
  %       r = floor(x ./ 2);
  %       y(InputDims{:}, b) = logical(x - 2*r);
  %       x = r;
  %     end
  
  %   if Transposed
  %     y = y';
  %   end
  
end

% Integer tests, in order of efficiency:
% x == fix(x); x == floor(x)
% mod(x,1) == 0
% rem(x,1) == 0
% % isequal(fix(x),x) % This only returns 1 value for the array.


% This was to use first singleton dimension for bits.  But more
% consistent to add extra dimension.
%   d = find(InputSize == 1, 1, 'first');
%   if ~isempty(d)
%     %     if d < nD
%     DimPerm = [1:d-1, nD, d+1:nD-1, d];
%     x = permute(x, DimPerm);
%   else
%     DimPerm = false;
%   end
% [...]
%   % Permute the dimensions back to original shape, with bits along
%   % ex-singleton dimension.
%   if any(DimPerm)
%     y = permute(y, DimPerm);
%   end;
