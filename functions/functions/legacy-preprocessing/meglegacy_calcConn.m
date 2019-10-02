function [ adjmat ] = meglegacy_calcConn( catdata, srate, varargin )
%MEGLEGACY_CALCCONN Calculate connectivity from a set of virtual sensor
%            time series, across trials.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  adjmat = meglegacy_calcConn(catdata, srate, ...)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%   adjmat -  adjacency / connectivity matrix describing connectivity over
%             time during the task between each sensor pair, for every time
%             point, for every frequency.
%             Has dimensions:
%               [sources] x [sources] x [time/samples] x [frequencies]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REQUIRED INPUTS:
%   catdata - concatenated virtual sensor timeseries, created in the virtual
%             sensor cells in meglegacy_1preprocess.m
%             Has dimensions:
%               [time/samples] x [trials] x [sensors]
%
%   srate -   sampling rate of the data. Can be determined from the .ds
%             header, or determined automatically in the virtual sensor
%             cells.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS:
%   metric -        (string) The metric or measure used to quantify the phase
%                   synchrony between the two signals.
%                   Valid inputs: 'pli', 'plv', 'wpli'
%
% Using this script to filter data
%   filterType -    ('fft', 'fir', or 'wavelet') Which time-frequency
%                   decomposition method to use to isolate specific
%                   frequencies for calculations.
%
% Parameters for different filter methods:
% %
% % Using wavelets:
% %   scales -        (vector) A vector of wavelet scales to use for computing the
% %                   wavelet transform. This is somewhat equivalent to
% %                   frequency bins. See `doc scal2freq` for more information.
% % 
% %   wname -         (string) Name of the wavelet function to use from the
% %                   MATLAB wavelet toolbox.
% % 
% % Using Fourier/FIR:
% %   freqs -         (matrix) List of frequency ranges to filter and perform
% %                   calculations on.
% %                   eg. [1,4; 4,7; 8,12;]
% % 
% % Using custom filter coefficients
% %   filtcoeffs -    (cell array) A cell array of filter coefficients to
% %                   convolve over data. The number of cells denotes the
% %                   number of frequency bands, and each cell should contain
% %                   a column vector of coefficients.
% % 
% %   filtfilt -      (1 or 0) Use (1) filtfilt two-pass function to filter
% %                   data or use (0) filter one-pass function.
% % 
% %   postFiltHilb -  (1 or 0) Perform Hilbert transform on data after
% %                   filtering. Set this to 0 if the filter already contains
% %                   a Hilbert transform, or set to 1 if not.
%
% Simeon Wong
% 2014 June 6

%% Process parameters
p = inputParser;

% Connectivity metric
addParamValue(p, 'metric', 'pli');

% Filter selection
addParamValue(p, 'filterType', 'fft');

% Wavelets
addParamValue(p, 'wname', 'cmor1-1');
addParamValue(p, 'scales', linspace(200,800,50));

% If providing frequency bins
addParamValue(p, 'freqs', []);

% If using your own filter coefficients
addParamValue(p, 'filtcoeffs', []);

addParamValue(p, 'filtfilt', 1);
addParamValue(p, 'postFiltHilb', 1);

parse(p, varargin{:});

%% Set things up

%%%%%
% Connectivity function
if strcmp(p.Results.metric, 'pli')
    % DOI: 10.1002/hbm.20346
    connfn = @(phasediff) abs(mean(sign(phasediff), 2));
elseif strcmp(p.Results.metric, 'plv')
    % DOI: 10.1002/(SICI)1097-0193(1999)8:4<194::AID-HBM4>3.0.CO;2-C
    connfn = @(phasediff) abs(mean(exp(1i .* phasediff), 2));
elseif strcmp(p.Results.metric, 'wpli')
    % DOI: 10.1016/j.neuroimage.2011.01.055
    % Equation 8
    connfn = @(phasediff) abs(sum(abs(phasediff) .* sign(phasediff), 2)) ./ sum(abs(phasediff), 2);
else
    error('%s is not a supported connectivity metric. (use pli, plv, or wpli)', p.Results.metric);
end

%%%%%
% Check for valid filter parameters
if isempty(p.Results.freqs) && isempty(p.Results.filtcoeffs)
    error('No filter parameters have been supplied!');
elseif ~isempty(p.Results.freqs) && ~isempty(p.Results.filtcoeffs)
    error('Only one of filter frequencies, or filter coefficients may be supplied!');
end

if p.Results.filtfilt == 1
    filterfn = @(b,x) filtfilt(b, 1, x);
    
    if p.Results.postFiltHilb == 0
        error('filtfilt cannot be used with filter coefficients combined with the Hilbert Transform. Consider setting parameter ''filtfilt'' to 0.');
    end
else
    filterfn = @(b,x) filter(b, 1, x);
end


%%%%%
% Get number of frequencies
if strcmp(p.Results.filterType, 'wavelet')
    nFreqs = length(p.Results.scales);
elseif ~isempty(p.Results.freqs)
    nFreqs = length(p.Results.freqs);
else
    nFreqs = length(p.Results.filtcoeffs);
end

[num_samples, num_trials, num_sensors] = size(catdata);

% Initialize result matrix
adjmat = zeros(num_sensors, num_sensors, num_samples, nFreqs);       % adjacency matrix over time

%% Calculations

if strcmp(p.Results.filterType, 'wavelet')
    %% Wavelets Calculations
    phasedata = zeros(nFreqs, num_samples, num_trials, num_sensors);
    parfor kk = 1:num_sensors
        for nn = 1:num_trials
            % get phase of signal using continuous wavelet transform
            phasedata(:,:,nn,kk) = angle(cwt(catdata(:,nn,kk), p.Results.scales, p.Results.wname));
        end
    end
    
    parfor f = 1:nFreqs
        paradjmat = zeros(num_sensors, num_sensors, num_samples);
        
        % take out required bit of phasedata
        % this becomes [samples] x [trials] x [sensors]
        parphasedata = squeeze(phasedata(f,:,:,:));
        
        % Loop through combinations of sensor pairs, and calculate connectivity
        for s1 = 1:num_sensors
            for s2 = s1+1:num_sensors
                %*** Calculate PLI ***
                % Get phase difference
                phasediff = parphasedata(:,:,s1) - parphasedata(:,:,s2);
                
                % Calculate connectivity, averaging across trials
                PLItimeseries = connfn(phasediff);
                
                paradjmat(s1,s2,:) = PLItimeseries;
                paradjmat(s2,s1,:) = PLItimeseries;    % assume bidirectional connectivity
            end
        end
        
        % fill result matrix
        adjmat(:,:,:,f) = paradjmat;
    end
    
else   
    %% Fourier-filter Calculations
    
    %%%%%
    % Loop over frequencies to perform calculations upon
    % Using parallel computing toolbox to speed up calculations
    parfor f = 1:nFreqs    %parfor
        
        % Create parallelized matrix to make MATLAB play nice
        paradjmat = zeros(num_sensors, num_sensors, num_samples);
        
        %*** Filter Data ***
        phasedata = zeros(num_samples, num_trials, num_sensors);
        
        % If frequency bins are supplied:
        if ~isempty(p.Results.freqs)
            
            % If filtering by FFT
            if strcmp(p.Results.filterType, 'fft')
                
                % Loop through virtual sensors, and trials and filter
                for kk = 1:num_sensors
                    for nn = 1:num_trials
                        % Filter and Hilbert Transform using fftFiltHilb
                        filtd = meglegacy_fftFiltHilb(catdata(:,nn,kk), ...
                            p.Results.freqs(f,1), p.Results.freqs(f,2), srate);
                        
                        % Extract phase
                        phasedata(:,nn,kk) = angle(filtd);
                    end
                end
                
                % If filtering by FIR least squares
            else
                error('not implemented');
            end
            
            % If filter coefficients are supplied
        else
            for kk = 1:num_sensors
                for nn = 1:num_trials
                    % Filter the data using filter function above
                    filtd = filterfn(p.Results.filtcoeffs{f}, catdata(:,nn,kk));
                    
                    % If Hilbert transform needs to be done for these
                    % coefficients, do it.
                    if p.Results.postFiltHilb == 1
                        filtd = hilbert(filtd);
                    end
                    
                    % Get phase from analytic (hilbert) signal
                    phasedata(:,nn,kk) = angle(filtd);
                end
            end
        end
        
        % Loop through combinations of sensor pairs, and calculate connectivity
        for s1 = 1:num_sensors
            for s2 = s1+1:num_sensors
                %*** Calculate PLI ***
                % Get phase difference
                phasediff = phasedata(:,:,s1) - phasedata(:,:,s2);
                
                % Calculate connectivity, averaging across trials
                PLItimeseries = connfn(phasediff);
                
                paradjmat(s1,s2,:) = PLItimeseries;
                paradjmat(s2,s1,:) = PLItimeseries;    % assume bidirectional connectivity
            end
        end
        
        adjmat(:,:,:,f) = paradjmat;
    end
end

end

