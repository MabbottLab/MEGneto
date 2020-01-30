function fcp_5_taskconnectivity(paths)

% FCP_5_TASKCONNECTIVITY estimates functional connectivity using phase-
% locking metrics (PLV, PLI) *or amplitude coupling* (to be added). 
% 
% NOTES:
%   - Will use the fcp_4 subject CSV to pull included participants.  
%   @JT: change AAL region number to automatic detection
%   @JT: add multi-condition handling
%
% INPUTS:
%   paths               =   struct defining paths to data, participant
%                           folders, analysis folders, config files, etc. 
%
% OUTPUTS:
%   adjmat:             participant adj matrix preserving results of each
%                       trial, saved into individual analysis folder
%   all_adjmat:         master matrix with all participants adj mats
%                       nodes x nodes x num_participants x num_freq_bands,
%                       saved into group analysis folder
%
% See also: 
%
% Last updated by: Julie Tseng, 2020-01-08
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.
%

%% SETUP

%%% PARTICIPANT IDS -------------------------------------------------------
% load config JSON with analysis parameters
config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp4';

% check for matched MRI and MEG data
subj_match  = ds_pid_match(paths,step);
ssSubjPath  = @(x) paths.(subj_match.pid{x});

%%% CONNECTIVITY FUNCTIONS ------------------------------------------------
if strcmp(config.connectivity.method,'pli')
    connfn     = @(H1,H2) abs(mean(sign(angle(H1)-angle(H2)), 1));
elseif strcmp(config.connectivity.method,'plv')
    connfn     = @(H1,H2) abs(mean(exp(1i .* (angle(H1)-angle(H2))), 1));
elseif strcmp(config.connectivity.method,'wpli')
    connfn     = @(H1,H2) abs(mean(abs(angle(H1)-angle(H2)).*(sign(angle(H1)-angle(H2))),1))/mean(angle(H1)-angle(H2));
elseif strcmp(config.connectivity.method,'coh')
    connfn     = @(H1,H2) abs(mean(H1.*conj(H2))./sqrt(mean(abs(H1).^2).*mean(abs(H2).^2)));
end

%%% OTHER USEFUL FUNCTIONS ------------------------------------------------
% Z-score
    prefilter = @(ts) zscore(ts);
% filter function
    filterfn = @(b,a,ts) filtfilt(b,a,ts);

%% CHECK DATA
% fileerror = false;

% %TO DO - change previous with following at some point  - SB
% % for cc = 1:length(p.condition)
%   for ss = 1:length(subj_match.ds)
%     %vspath = p.paths.vs(char(p.subj.ID{ss}),
%     %p.beamformer.bf.method,p.condition{cc});  to change later
%     vspath = [ssSubjPath(ss) '/AAL_beamforming_results'];
%     if ~exist(vspath, 'file')
%       fileerror = true;
%       disp([vspath ' doesn''t exist!']);
%     else
%       vsdata = whos('-file', vspath, 'catmatrix', 'srate');
%       
%       if length(vsdata) < 2
%         fileerror = true;
%         disp([vspath ' is corrupted!']);
%       end
%     end
%   end
% % end

% if fileerror
%   error('Some source files seem to have issues!');
% end

%% MAKE FILTERS

%%% LOAD SAMPLE PPT -------------------------------------------------------
load([ssSubjPath(1) '/AAL_beamforming_results'], '-mat') %SB - add cond

%%% MAKE FILTER -----------------------------------------------------------
maxn = 0;
for fq = 1:length(config.connectivity.filt_freqs)

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

% setup band names and master adjacency matrix
band_names = ["theta", "alpha", "beta", "lowgamma", "highgamma"];
all_adjmat = nan(90, 90, length(subj_match.ds), length(config.connectivity.filt_freqs));

%for cc = 1:length(p.condition)
  % fprintf('Running condition %s... \n', p.condition{cc});
  
    for ss = 1:length(subj_match.ds)
        fprintf('\n\n==================================\nSUBJECT: %s\n', subj_match.pid{ss});

        fprintf('Loading participant %s... \n',subj_match.pid{ss});

%%% LOAD VIRTUAL SENSOR DATA ----------------------------------------------
        load([ssSubjPath(ss) '/AAL_beamforming_results'], '-mat'); 

        % define some dimensions
        num_samples = size(catmatrix, 1);
        num_trials  = size(catmatrix, 2);
        num_sources = size(catmatrix, 3);

%%% INITIALIZE PARTICIPANT ADJACENCY MATRIX -------------------------------
%   Dimensions = [sources] x [sources] x [trials] x [freq. band]
        adjmat = nan(num_sources, num_sources, num_trials,length(config.connectivity.filt_freqs)); % in resting state

%%% RUN CONNECTIVITY ANALYSIS ---------------------------------------------
    %%% FOR EACH FREQUENCY BAND
        for fq = 1:length(config.connectivity.filt_freqs) 
          fprintf('Analyzing the %s band!\n', band_names(fq)); 

          % initialize the temp adjmat
          p_adjmat = nan(num_sources, num_sources, num_trials); % same as RestingSate     
          H_data = complex(zeros(num_samples, num_trials, num_sources)); 
          
          %%% FOR EACH TRIAL
          fprintf('Processing trial ');
          for tt = 1:num_trials
            fprintf('\b\b\b\b%d...', tt);
            %%% FOR EACH AAL NODE/SOURCE
            for kk = 1:num_sources
              % mean center or z-score the timeseries before filtering
              ts = prefilter(catmatrix(:,tt,kk));          
              % filter data, calculate hilbert transform, get instantaneous phase
              if (max(length(fir_coef{fq})-1,length(1)-1)*3) < length(ts)          
                H_data(:,tt,kk) = hilbert(filtfilt(fir_coef{fq}, 1, ts));
              else
                H_data(:,tt,kk) = hilbert(filter(fir_coef{fq}, 1, ts));
              end
            end
          end
          fprintf('\n');
%%% CALCULATE CONNECTIVITY ------------------------------------------------
          fprintf('Onto the connectivity calculations!\n')
          for aa = 1:num_sources
            for bb = aa+1:num_sources
              p_adjmat(aa,bb,:) = connfn(H_data(:,:,aa), H_data(:,:,bb));
              p_adjmat(bb,aa,:) = p_adjmat(aa,bb,:);
            end
          end
          fprintf('Done this band.')
%%% COPY TEMP MATRIX INTO SUBJECT-MATRIX ----------------------------------
          adjmat(:,:,:,fq) = p_adjmat;
        end % repeat for each frequency band

%%% SAVE PARTICIPANT'S ADJACENCY MATRIX -----------------------------------
        fprintf('Saving the mat file...');
        save([ssSubjPath(ss) '/fcp_5_adjmat_' config.connectivity.method '.mat'],'adjmat','-mat','-v7.3')

%%% CALCULATE AND RECORD AVG ACROSS TRIALS FOR THIS PARTICIPANT
        all_adjmat(:,:,ss,:) = squeeze(mean(adjmat, 3));
    end
    
%%% SAVE MASTER ADJACENCY MATRIX WITH ALL PARTICIPANTS
save([paths.anout_grp '/fcp_5_allParticipants_adjmats_' config.connectivity.method '.mat'],'all_adjmat','-v7.3');

fprintf('\n\n============== Finished Processing ====================\n');  
