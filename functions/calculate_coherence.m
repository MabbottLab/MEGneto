function out = calculate_coherence(x, y, fs, T, foi)

% x, y:     time x trial matrix, e.g., 2400 x 98
% dt:       sampling interval in seconds
% T:        total duration of data
% df:       1/T; frequency resolution
% fNQ:      1/dt/2; nyquist frequency
% foi:      frequency band

df = 1/T;
dt = 1/fs;
fNQ = fs/2;

num_trials = size(x, 2);

for tr = 1:num_trials
    Sxx(tr,:) = 2*dt^2/T*(fft(x(:,tr)').*conj(fft(x(:,tr)')));
    Syy(tr,:) = 2*dt^2/T*(fft(y(:,tr)').*conj(fft(y(:,tr)')));
    Sxy(tr,:) = 2*dt^2/T*(fft(x(:,tr)').*conj(fft(y(:,tr)')));
end

Sxx = nanmean(Sxx(:,1:(fs/2+1)),1);
Syy = nanmean(Syy(:,1:(fs/2+1)),1);
Sxy = nanmean(Sxy(:,1:(fs/2+1)),1);

coh = abs(Sxy)./(sqrt(Sxx).*sqrt(Syy));

faxis = 0:df:fNQ;
include_these = (faxis >= foi(1) & faxis <= foi(2));
out = nanmean(coh(include_these));

end