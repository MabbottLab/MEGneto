function Step4a_FrequencyAnalysis(config, pid)

% docs

%% load setup

this_output = [config.meta.project_path '/' config.meta.analysis_name '/' pid];
load([this_output '/out_struct.mat'])

%% RUN FREQUENCY ANALYSIS

    % load data
    load([this_output '/step3_data_roi.mat'])
        
    %%% Timewindow analysis
    if config.step4a.type == "power"
        % select data from each timewindow of interest
        cfg = [];
        cfg.latency = config.step4a.toi;
        data_toi        = ft_selectdata(cfg, data_roi);
        
        % select data from baseline
        cfg.latency = config.step4a.boi;
        data_boi        = ft_selectdata(cfg, data_roi);
        
        % generate power spectra for toi and boi
        cfg             = [];
        cfg.output      = 'pow';
        cfg.taper       = 'hanning';
        cfg.foi         = config.step4a.foi; % frequencies of interest
        cfg.method      = 'mtmfft';
        freq_toi        = ft_freqanalysis(cfg, data_toi);
        freq_boi        = ft_freqanalysis(cfg, data_boi);
        
        % instantiate resulting struct and retain helpful info
        freq = keepfields(freq_toi_avg, {'label', 'cfg', 'dimord'});
        freq.latency = [freq_toi.cfg.previous.latency; ...
                        freq_boi.cfg.previous.latency];
        freq.freq = [];
        freq.powspctrm = [];
        cfg = [];
        for fband = 1:length(config.step4a.freqbands)
            cfg.frequency   = config.step4a.freqbands(fband,:);
            cfg.avgoverfreq = 'yes';
            cfg.nanmean     = 'yes';
            freq_toi_avg    = ft_selectdata(cfg, freq_toi);
            freq_boi_avg    = ft_selectdata(cfg, freq_boi);
            
            freq.freq       = [freq.freq; ...
                               freq_toi_avg.freq freq_boi_avg.freq];
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
        % actual timewindow analysis
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
        freq = keepfields(freq_bl, {'label', 'dimord', 'time', 'cfg'});
        freq.freq = [];
        freq.powspctrm = [];
        for fband = 1:length(config.step4a.freqbands)
            cfg = [];
            cfg.avgoverfreq = 'yes';
            cfg.frequency = config.step4a.freqbands(fband,:);
            cfg.nanmean = 'yes';
            this_avg = ft_selectdata(cfg, freq_bl);
            
            freq.freq = [freq.freq; this_avg.freq];
            freq.powspctrm = cat(2, freq.powspctrm, this_avg.powspctrm);
        end
        
        save([this_output '/step4a_freq_tfr.mat'], 'freq')
    end 
    
end