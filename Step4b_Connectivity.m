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
    
    if ~isfield(config.step4b, 'conditions')
        config.step4b.conditions = {"all", 1:length(data_roi.trials)};
    else
        config.step4b.conditions(:,2) = ...
            arrayfun(@(x) find(data_roi.trialinfo == x), ...
            [config.step4b.conditions{:,2}], 'UniformOutput', false)';
    end
    
    conn.conditions = config.step4b.conditions;
    
    for fq = 1:length(config.step4b.freqbands) 
        cfg = [];
        cfg.bpfilter = 'yes';
        cfg.bpfreq = config.step4b.freqbands(fq,:);
        data_roi_filt = ft_preprocessing(cfg, data_roi);

        for cond = 1:size(config.step4b.conditions, 1)
              % select relevant time period
              cfg = [];
              cfg.latency = config.step4b.toi;
              cfg.trials = config.step4b.conditions{cond,2};
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
                    conn.(sprintf('%sspctrm',cfg.method))(cond,:,:,fq) = ...
                        squeeze(max(this_conn.(sprintf('%sspctrm', cfg.method)),[], 3));
                else % otherwise, use mean
                    conn.(sprintf('%sspctrm',cfg.method))(cond,:,:,fq) = ...
                        squeeze(mean(this_conn.(sprintf('%sspctrm', cfg.method)), 3));
                end
        end
    end

    save([this_output '/step4b_conn.mat'], 'conn');

end
