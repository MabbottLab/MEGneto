function permutations(paths)

% Permutations estimates functional connectivity (i.e., analyzes 
% the synchrony of signals from two regions). 
% Currently supports the following connectivity metrics: 
%       - 'plv' (phase locking value)
%       - 'pli' (phase lag index)
%       - 'wpli' (weighted phase lag index)
%       - 'wpli_debiased' (as above, but debiased)
%       - 'coh' (coherence)
% The strings listed above should be indicated in the config JSON file
% exactly as presented. 
%
% INPUTS:
%   paths               =   struct defining paths to data, participant
%                           folders, analysis folders, config files, etc. 
%
% OUTPUTS:
%   conn_mat:             participant connectivity matrix preserving 
%                         results of each trial, saved into individual 
%                         analysis folder
%   all_conn_mat:         master matrix with all participants conn mats
%                         nodes x nodes x num_participants x num_freq_bands,
%                         saved into group analysis folder
%
% See also: 
%
% Last updated by: Julie Tseng, 2020-01-08
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.
%
%% SET UP LOGGING FILE

right_now = clock;
log_filename = [paths.conf_dir '/log_' sprintf('%02.f%02.f%02.f', right_now(1:3))];
diary(log_filename)

fprintf('\n\n%02.f:%02.f:%02.f       Now running **%s**.\n', ...
    right_now(4:6), mfilename)

%% SETUP

%%% PARTICIPANT IDS -------------------------------------------------------
% load config JSON with analysis parameters
config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp5';

% check is user indicated that connectivity results should be generated
if ~config.connectivity.include
    error('Your input to the JSON config file indicates you do not wish to run fcp_5_taskconnectivity. If you wish to run this function, please change the config.connectivity.include field to 1.')
end

% check for matched MRI and MEG data
subj_match  = ds_pid_match(paths,step);
ssSubjPath  = @(x) paths.(subj_match.pid{x});

%%% OTHER USEFUL FUNCTIONS ------------------------------------------------
% Z-score
    prefilter = @(ts) zscore(ts);
% filter function
    filterfn = @(b,a,ts) filtfilt(b,a,ts);

%% MAKE FILTERS

%%% LOAD SAMPLE PPT -------------------------------------------------------
% used to get the sample rate
% used compiled CAT matrices here
try
    load([ssSubjPath(1) '/compiledCATmatrix.mat'], '-mat');
catch
    load([ssSubjPath(1) '/AAL_beamforming_results.mat'], '-mat');
end

%%% MAKE FILTER -----------------------------------------------------------
maxn = 0;
for fq = 1:length(config.connectivity.filt_freqs) % for each frequency band specified in the JSON config file

    % calculate filter coefficients
    filtband = config.connectivity.filt_freqs(fq,:);
  
    % estimate filter order
    transition = mean(filtband) * 0.2;
    freqs = [filtband(1) - transition, filtband(1), filtband(2), filtband(2) + transition];
    fir_n = kaiserord(freqs, [0 1 0], [0.1 0.05 0.1], srate);
    maxn = max(maxn, fir_n);
  
    % make filter coefficients
    fir_coef{fq} = fir1(fir_n, filtband*2/srate, 'bandpass');

%   if (max(length(fir_coef{fq})-1,length(1)-1)*3) > length(ts)
%       fprintf('Warning: For fiter [%d %d] - zero phase filtering can not be used.\n' ,p.filt_freqs(fq,:)) 
%   end
  
end

fprintf('Max filter length: %d samples = %.4f sec.\n', maxn, maxn/srate);

%% RUN CONNECTIVITY ANALYSIS

% setup band names and master connectivity matrix
band_names = config.connectivity.freq_names;

% get number of ROIs or supra-ROIs
ROIs = config.connectivity.ROIs;
if ~isempty(ROIs{1})
    num_sources = length(ROIs); % if collapsed
else
    num_sources = size(catmatrix, 3); % if not
end

% tmp_conn_mat = nan(num_sources, num_sources, 1, 5);
all_conn_mat = nan(num_sources, num_sources, ... % num_nodes x num_nodes x ...
                   500, ...                      % num permutations ...
                   length(config.connectivity.filt_freqs));    % num_freq_bands

%for cc = 1:length(p.condition)
  % fprintf('Running condition %s... \n', p.condition{cc});
    for ss = 1:12
        right_now = clock;
        fprintf('%02.f:%02.f:%02.f       Working on participant %s!\n', ...
            right_now(4:6), num2str(ss))

%%% LOAD VIRTUAL SENSOR DATA ----------------------------------------------
        try
            catmatrix = load([feval(ssSubjPath, ss) '/compiledCATmatrix'], '-mat', 'catmatrix'); 
            catmatrix = catmatrix.catmatrix;
        catch
            catmatrix = load([feval(ssSubjPath, ss) '/AAL_beamforming_results.mat'], '-mat', 'catmatrix');
            catmatrix = catmatrix.catmatrix;
        end
    
        % define some dimensions
        num_samples = size(catmatrix, 1);
        num_trials  = size(catmatrix, 2);
        
        % collapse across ROIs if indicated
        if ~isempty(ROIs{1})
            catmatrix_collapsed = cellfun(@(x) nanmean(catmatrix(:,:,x), 3), ROIs, 'UniformOutput', false);
            catmatrix = cat(3, catmatrix_collapsed{:}); clear catmatrix_collapsed;
        end

%%% INITIALIZE PARTICIPANT CONNECTIVITY MATRIX -------------------------------
%   Dimensions = [sources] x [sources] x [freq. band]
%         conn_mat  = nan(num_sources, num_sources, length(config.connectivity.filt_freqs)); 
        conn_mat  = nan(num_sources, num_sources, 1); 
        data    = [];
        time_info = config.task.trialdef.parameters.tEpoch;
        for src = 1:num_sources
            data.label{src} = sprintf('ROI%d', src);
        end

%%% RUN CONNECTIVITY ANALYSIS ---------------------------------------------
    %%% FOR EACH FREQUENCY BAND
    
%           for fq = 1:length(config.connectivity.filt_freqs) % for loop for frequencies
            fq = 4;
          fprintf('Analyzing the %s band!\n', band_names{fq}); 
          counter = 1;
          %%% FOR EACH TRIAL
          for tt = 1:num_trials
              fprintf('Processing trial %d...\n', tt);
              %           tmp = sum(~isnan(catmatrix(:, tt, 1)));
              tmp = catmatrix(:, tt, 1);
              remove_indx = isnan(tmp);
              tmp(remove_indx) = [];
              data.time{counter} = (0.400 + (1/srate)):(1/srate):(0.4+length(tmp)/srate);
              %%% FOR EACH OF THE NODES/SOURCES
              for kk = 1:num_sources
                  % mean center or z-score the timeseries before filtering
                  tmp = catmatrix(:, tt, kk);
                  remove_indx = isnan(tmp);
                  tmp(remove_indx) = [];
                  ts = prefilter(tmp); % prepared specified portion of catmatrix
                  % filter data, calculate hilbert transform, get instantaneous phase
                  if (max(length(fir_coef{fq})-1,length(1)-1)*3) < length(ts)
                      data.trial{counter}(kk,:) = filtfilt(fir_coef{fq}, 1, ts);
                  else
                      data.trial{counter}(kk,:) = filter(fir_coef{fq}, 1, ts);
                  end
              end
              counter = counter+1;
          end
         
%%% CALCULATE CONNECTIVITY ------------------------------------------------
            fprintf('Onto the connectivity calculations!\n')
            %% do 5000 permutations of random stuff, then see if the social content
            % is significant. save each run of permutation, then plot by
            % pairwise ROI, see where social content lies.
            % don't forget to make if and elseif statements for ppts 8, 16,
            % and 17
%             load([ssSubjPath(ss) '/socialcontent.mat']);
%             social          = find(social);
            cfg             = []; % set up config for connectivity calculation
            cfg.method      = 'mtmfft';
            cfg.output      = 'powandcsd';
            cfg.channel     = 'all';
            cfg.trials      = 'all';
            cfg.keeptrials  = 'yes';
            cfg.taper       = 'hanning';
            cfg.foilim      = config.connectivity.filt_freqs(fq,:); % calculate connectvity'
            cfg.pad         = 5;
            % GET THE CROSS SPECTRUM
            freq            = ft_freqanalysis(cfg, data); 

            parfor jj = 1:1000
                conn_mat = zeros(100, 100, 1);
                content = randperm(num_trials, 52);
                % COMPUTE CONNECTIVITY METRIC
                if strcmp(config.connectivity.method, 'pli')
                    % need to trick fieldtrip if we're using PLI, as they don't
                    % have it implemented in their ft_connectivityanalysis
                    conn.labelcmb      = freq.labelcmb;
                    conn.dimord        = 'chancmb_freq';
                    conn.freq          = freq.freq;
                    conn.cfg           = freq.cfg;
                    conn.plispctrm     = squeeze(abs(mean(sign(freq.crsspctrm), 1)));
                else
                    cfg             = []; % set up config for computing connectivity metric
                    cfg.method      = config.connectivity.method;
                    cfg.trials      = content;
                    conn            = ft_connectivityanalysis(cfg, freq); % compute connectivity metric
                end
                
                % RESHAPE INTO SOURCE X SOURCE CONN MAT AND STORE
                conn            = ft_checkdata(conn, 'cmbrepresentation', 'full');
                if ~strcmp(config.connectivity.method, config.connectivity.collapse_band) % default to max across band
                    conn_mat(:,:,1) = squeeze(max(conn.(sprintf('%sspctrm', config.connectivity.method)),[], 3));
                else % otherwise, use mean
                    conn_mat(:,:,1) = squeeze(mean(conn.(sprintf('%sspctrm', config.connectivity.method)), 3));
                end
                
                fprintf('Done this band.')
%                 all_conn_mat(:,:,jj) = conn_mat;
%             end
              ppt(:,:,jj,ss) = conn_mat;
            end
    tmp_conn_mat = ppt;
%     end % repeat for each frequency band
    
    %%% SAVE PARTICIPANT'S CONNECTIVITY MATRIX -----------------------------------
%     fprintf('Saving the mat file...\n');
%     save([ssSubjPath(ss) '/fcp_5_controlconn_mat_' config.connectivity.method '.mat'],'conn_mat','-mat','-v7.3')
    
    %%% CALCULATE AND RECORD AVG ACROSS TRIALS FOR THIS PARTICIPANT
    
% end

%%% SAVE MASTER CONNECTIVITY MATRIX WITH ALL PARTICIPANTS
save([paths.anout_grp '/fcp_5_lowgamma_raw_control_mats_' config.connectivity.method '.mat'],'tmp_conn_mat','-v7.3');

%% turn off diary
right_now = clock;
fprintf('%02.f:%02.f:%02.f ============== Finished Processing ====================\n', ...
    right_now(4:6))
diary off

%% let the users know that connectivity analysis is complete
sendEmail("connectivity", string(config.contact));
end