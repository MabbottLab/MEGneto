function [WPLIMatrix] = pairwise_wPLI_wavelet(x)
% Input: wavelet coefficients [time_points channels]

% x = hilbert(data);

[rows,cols] = size(x);
if rows < cols,
    x = x';
end

H_imag=imag(x);
H_real=real(x); 

nchan = size(H_imag,2);
WPLIMatrix = zeros(nchan,nchan);

for col=1:size(H_imag,2)
    imagX = bsxfun(@times,H_imag,H_real(:,col))-bsxfun(@times,H_real,H_imag(:,col));
    expectedImag        = mean(imagX);
    expectedAbs         = mean(abs(imagX));
    WPLIMatrix(col,:)  = abs(expectedImag)./expectedAbs;  
end

WPLIMatrix(1:nchan+1:end) = 1; 

return
