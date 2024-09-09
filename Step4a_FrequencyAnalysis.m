function Step4a_FrequencyAnalysis(config, pid, visitnum)

% Step4a_FrequencyAnalysis uses spectral analysis on time-frequency 
% representations of data to test hypotheses based on spectral power. 
% The virtual sensor data from the beamforming step is loaded in and 
% frequency analysis is performed on sliding timewindows of the data. Thus, 
% for each subject and each interpolated atlas region, a power spectrum is 
% calculated and corrected to a baseline to control for general/random spikes in power. 
% This can be configured to be an overall power analysis, or a TFR (sliding 
% time window) power analysis. 
%
% INPUTS:
%   > config: 
%       struct, configured in your "main" script with
%       all analysis parameters and options and paths
%   > pid:
%       string, participant ID used to build the paths
%       to relevant data files and output folders
%   > ds_path:
%       string, full path to MEG *.ds folder for this participant
%   > visitnum:
%       int, visit number if longitudinal, optional arg)
%
% OUTPUTS:
%   > step1_data_clean.mat:
%       fieldtrip style data object with data epoched into trials and 
%       noisy trials (head motion, artifact) identified + rejected
%   > out_struct.mat:
%       a MATLAB structure with meta-information about the pipeline step
%       run, e.g.: number of trials epoched, left after rejection
%   > plot_markers.png:
%       a PNG image visualizing the markers found in the *.ds file, that
%       can be used as a quick inspection of task markers
%   > headmotion_[date].png:
%       a PNG image outputted by the HeadMotionTool visualizing head motion
%       across the timeseries

% Last updated by: Julie Tseng, 2024-09-09
%   This file is part of MEGneto, see https://github.com/MabbottLab/MEGneto
%   for the documentation and details.

%% SETUP

if exist('visitnum', 'var')
    this_output = [config.meta.project_path '/' config.meta.analysis_name '/' ...
        pid '/' sprintf('ses-%.2d', visitnum)]; % indicate subject-specific output folder path
else
    this_output = [config.meta.project_path '/' config.meta.analysis_name '/' pid]; 
end

load([this_output '/out_struct.mat'])

%% RUN FREQUENCY ANALYSIS

    % load data
    load([this_output '/step3_data_roi.mat'])
        
    %%% Timewindow analysis
    if config.step4a.type == "power"
        % select data from each timewindow of interest
        cfg             = [];
        cfg.latency     = config.step4a.toi;
        data_toi        = ft_selectdata(cfg, data_roi);
        
        % select data from baseline
        cfg.latency     = config.step4a.boi;
        data_boi        = ft_selectdata(cfg, data_roi);
        
        % generate power spectra for toi and boi
        cfg             = [];
        cfg.output      = 'pow';
        cfg.taper       = 'hanning';
        cfg.foi         = config.step4a.foi; % frequencies of interest
        cfg.method      = 'mtmfft';
        freq_toi        = ft_freqanalysis(cfg, data_toi); % timewindow of interest
        freq_boi        = ft_freqanalysis(cfg, data_boi); % baseline
        
        % instantiate resulting struct and retain helpful info
        freq = keepfields(freq_toi_avg, {'label', 'cfg', 'dimord'});
        freq.latency    = [freq_toi.cfg.previous.latency; ...
                            freq_boi.cfg.previous.latency];
        freq.freq       = [];
        freq.powspctrm  = [];
        cfg             = [];

        % for each frequency band
        for fband = 1:length(config.step4a.freqbands)
            % set the frequencies
            cfg.frequency   = config.step4a.freqbands(fband,:);
            cfg.avgoverfreq = 'yes';
            cfg.nanmean     = 'yes';
            freq_toi_avg    = ft_selectdata(cfg, freq_toi);
            freq_boi_avg    = ft_selectdata(cfg, freq_boi);
            
            % append the frequencies of the toi/boi to the resulting freq
            % structure for bookkeeping
            freq.freq       = [freq.freq; ...
                               freq_toi_avg.freq freq_boi_avg.freq];
                           
            % choose your particular baseline type
            if strcmp(config.step4a.baselinetype, 'absolute')
                freq.powspctrm  = [freq.powspctrm ...
                                   freq_toi_avg.powspctrm - freq_boi_avg.powspctrm];
            elseif strcmp(config.step4a.baselinetype, 'relative')
                freq.powspctrm = [freq.powspctrm ...
                                   freq_toi_avg.powspctrm ./ freq_boi_avg.powspctrm];
            elseif strcmp(config.step4a.baselinetype, 'relchange')
                freq.powspctrm = [freq.powspctrm ...
                                   (freq_toi_avg.powspctrm - freq_boi_avg.powspctrm) ./ freq_boi_avg.powspctrm];
            end
        end
        save([this_output '/step4a_freq_overallpower.mat'], 'freq')
        
    elseif config.step4a.type == "tfr"
        % timewindow analysis
        cfg             = [];
        cfg.method      = 'mtmconvol';
        cfg.taper       = 'hanning';
        cfg.foi         = config.step4a.foi;
        cfg.t_ftimwin   = 4 ./cfg.foi;
        cfg.toi         = config.step4a.toi; % time of interest
        freq_bl         = ft_freqanalysis(cfg, data_roi); % perform frequency analysis
        
        % baselining
        cfg             = [];
        cfg.baseline    = config.step4a.boi;
        cfg.baselinetype = config.step4a.baselinetype;
        freq_bl         = ft_freqbaseline(cfg, freq_bl);
    
        % avg by frequency band
        freq        = keepfields(freq_bl, {'label', 'dimord', 'time', 'cfg'});
        freq.freq   = [];
        freq.powspctrm = [];
        
        % for each frequency band
        for fband = 1:length(config.step4a.freqbands)
            cfg                 = [];
            cfg.avgoverfreq     = 'yes';
            cfg.frequency       = config.step4a.freqbands(fband,:);
            cfg.nanmean         = 'yes';
            this_avg            = ft_selectdata(cfg, freq_bl);

            % bookkeeping
            freq.freq = [freq.freq; this_avg.freq];
            
            % resulting powspctrm is a multidimensional matrix with one of
            % the dimensions being frequency band
            freq.powspctrm = cat(2, freq.powspctrm, this_avg.powspctrm);
        end
        save([this_output '/step4a_freq_tfr.mat'], 'freq')
    end 
    
end