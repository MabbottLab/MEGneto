function Step4b_Connectivity(config, pid, visitnum)
%
% docs
%
%% SETUP

if exist('visitnum', 'var')
    this_output = [config.meta.project_path '/' config.meta.analysis_name '/' pid '/' sprintf('ses-%.2d', visitnum)]; % indicate subject-specific output folder path
else
    this_output = [config.meta.project_path '/' config.meta.analysis_name '/' pid]; 
end

load([this_output '/out_struct.mat'])
load([this_output '/step3_data_roi.mat'])

%% RUN CONNECTIVITY ANALYSIS

%%% RUN CONNECTIVITY ANALYSIS ---------------------------------------------
    %%% FOR EACH FREQUENCY BAND
    conn = [];
    conn.dimord = 'chan_chan_freq';
    conn.label = data_roi.label;
    
    for fq = 1:length(config.step4b.freqbands) 
          cfg = [];
          cfg.bpfilter = 'yes';
          cfg.bpfreq = config.step4b.freqbands(fq,:);
          data_roi_filt = ft_preprocessing(cfg, data_roi);
          
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
            this_conn            = ft_connectivityanalysis(cfg, freq_filt); % compute connectivity metric
            
            % RESHAPE INTO SOURCE X SOURCE CONN MAT AND STORE
            this_conn            = ft_checkdata(this_conn, 'cmbstyle', 'full');
            conn.freq(fq) = mean(this_conn.freq);
            
            % collapse frequencies within band
            if ~strcmp(config.step4b.bandavgmethod, 'max') % default to max across band
                conn.(sprintf('%sspctrm',cfg.method))(:,:,fq) = ...
                    squeeze(max(this_conn.(sprintf('%sspctrm', cfg.method)),[], 3));
            else % otherwise, use mean
                conn.(sprintf('%sspctrm',cfg.method))(:,:,fq) = ...
                    squeeze(mean(this_conn.(sprintf('%sspctrm', cfg.method)), 3));
            end
    end
    
    save([this_output '/step4b_conn.mat'], 'conn');

end
