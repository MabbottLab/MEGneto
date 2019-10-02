function [ dataOut ] = meglegacy_fftFiltHilb( dataIn, f1, f2, srate )
%MEGLEGACY_FFTFILTHILB Filters data using an FFT and calculates hilbert transform
%
% Uses an FFT and inverse FFT to filter and calculate hilbert transform of
% the input data set in one step.
% A bandpass filter is applied to the data between the specified
% frequencies.
%
% Arguments:
%   dataIn- Input data. Can be column or row vector, and must be consist of
%           only a single trial
%   f1-     Bottom end of bandpass frequency band in Hz
%   f2-     Top end of bandpass frequency band in Hz
%   srate-  Sampling rate of data in Hz
% 
% References:
%   - eegfiltfft    (simplified version of this function with added Hilbert)
% 
% Simeon Wong
% 31 May 2012
% Originally intended for preterm pipeline (Sam Doesburg)
%
% Modified 17 May 2013 for epilepsy language expression pipeline

    dataSize = length(dataIn);
    
    freqs = 0:srate/(dataSize-2):srate/2;
    
    lowIndices = find((freqs < f1));
    highIndices = find((freqs > f2));
   
   
    dfft = fft(dataIn);
    
    dfft([lowIndices highIndices highIndices(end):end]) = 0;
     
    dfft = dfft .* 2;
   
    dataOut = ifft(dfft);

end

