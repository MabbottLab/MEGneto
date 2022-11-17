function Step4b_Connectivity(config, pid)
%
% docs
%
%% SETUP

this_output = [config.meta.project_path '/' config.meta.analysis_name '/' pid];
load([this_output '/out_struct.mat'])
load([this_output '/step3_data_roi.mat'])

%% MAKE FILTERS

for fq = 1:length(config.step4b.freqbands) % for each frequency band specified in the JSON config file
    filtband = config.step4b.freqbands(fq,:);
    transition = mean(filtband) * 0.2;
    
    bpfilt{fq} = designfilt('bandpassfir', ...
                                'StopbandFrequency1', filtband(1) - transition, ...
                                'PassbandFrequency1', filtband(1), ...
                                'PassbandFrequency2', filtband(2), ...
                                'StopbandFrequency2', filtband(2) + transition, ...
                                'DesignMethod', 'kaiserwin', ...
                                'PassbandRipple', 0.5, ...
                                'StopbandAttenuation1', 35, ...
                                'StopbandAttenuation2', 35, ...
                                'SampleRate', data_roi.fsample);
end

%% RUN CONNECTIVITY ANALYSIS

%%% RUN CONNECTIVITY ANALYSIS ---------------------------------------------
    %%% FOR EACH FREQUENCY BAND
    freq = [];
    freq.dimord = 'chan_chan_freq';
    freq.label = data_roi.label;
    
    for fq = 1:length(bpfilt) 
          data_roi_filt = keepfields(data_roi, {'time', 'fsample', 'label'});
          for tt = 1:length(data_roi.trial)
              try
                this_filtered = arrayfun(@(x) filtfilt(bpfilt{fq}, ...
                                                        data_roi.trial{tt}(x,:)), ...
                                                        1:length(data_roi.label), ...
                                                        'UniformOutput', false);
              catch
                  this_filtered = arrayfun(@(x) filter(bpfilt{fq}, ...
                                                        data_roi.trial{tt}(x,:)), ...
                                                        1:length(data_roi.label), ...
                                                        'UniformOutput', false);
              end
              data_roi_filt.trial{tt} = cat(1,this_filtered{:});
          end
                 
          
          % select relevant time period
          cfg = [];
          cfg.latency = config.step4b.toi;
          data_roi_filt_toi = ft_selectdata(cfg, data_roi_filt);
          
%%% CALCULATE CONNECTIVITY ------------------------------------------------
            fprintf('Onto the connectivity calculations!\n')
            cfg             = []; % set up config for connectivity calculation
            cfg.method      = 'mtmfft';
            cfg.output      = 'powandcsd';
            cfg.channel     = 'all';
            cfg.trials      = 'all';
            cfg.keeptrials  = 'yes';
            cfg.taper       = 'hanning';
            cfg.foilim      = config.step4b.freqbands(fq,:);
            freq_filt       = ft_freqanalysis(cfg, data_roi_filt_toi); 

            % COMPUTE CONNECTIVITY METRIC
            cfg             = []; % set up config for computing connectivity metric
            cfg.method      = config.step4b.connmethod;
            conn            = ft_connectivityanalysis(cfg, freq_filt); % compute connectivity metric
            
            % RESHAPE INTO SOURCE X SOURCE CONN MAT AND STORE
            conn            = ft_checkdata(conn, 'cmbstyle', 'full');
            freq.freq(fq) = mean(conn.freq);
            
            % collapse frequencies within band
            if ~strcmp(config.step4b.bandavgmethod, 'max') % default to max across band
                freq.(sprintf('%sspctrm',cfg.method))(:,:,fq) = ...
                    squeeze(max(conn.(sprintf('%sspctrm', cfg.method)),[], 3));
            else % otherwise, use mean
                freq.(sprintf('%sspctrm',cfg.method))(:,:,fq) = ...
                    squeeze(mean(conn.(sprintf('%sspctrm', cfg.method)), 3));
            end
    end
    

end
