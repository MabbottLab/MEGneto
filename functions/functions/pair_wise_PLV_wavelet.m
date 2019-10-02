function [PLV,meanphi] = pair_wise_PLV_wavelet(data)
%function PLV = pairwisePLV(data)
% Compute Phase Locking Value (PLV) for all pairs of channels
% Data: timepoints x channels
%
%
% Example:
% x = sin(2*10*pi*t);
% y = sin(2*10*pi*t);
% x = x + randn(size(x))*0.1;
% y = y + randn(size(y))*0.1;
% [PLV,meanphi] = pair_wise_PLV_wavelet([x' x'+i*y'])
% Output:
%   PLV =
%   1.0000    0.4825
%    0.4825    1.0000
%   meanphi =
%     0   -90
%    90     0
% Vasily - May 06,2011

[rows,cols] = size(data);
if rows < cols,
    data = data';
end

cols = size(data,2);

% initialize matrix for efficient memory use
PLV = eye(cols);
meanphi = zeros(cols);

for ii = 1:cols-1
  for jj = (ii+1):cols

    x = data(:,ii);
    y = data(:,jj);

    % take the difference in instantaneous phase between the signals
    delta = angle(x.*conj(y));
    mean_delta = meanangle(degrees(delta),1);
    
    % construct a series of 2-D vectors on unit circle, one for each datapoint
    vectors = [cos(delta) sin(delta)];
    mv = mean(vectors,1);
    
    %phi = atan2(mv(2),mv(1));
    % the PLI is the length of the mean of those vectors
    PLV(ii,jj) = norm(mv);    
    PLV(jj,ii) = PLV(ii,jj);
    meanphi(ii,jj) = mean_delta;
    meanphi(jj,ii) = -mean_delta;
    
  end
end

   % keyboard
    
    return