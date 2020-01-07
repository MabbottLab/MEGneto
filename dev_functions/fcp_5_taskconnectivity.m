function fcp_5_taskconnectivity(paths)
% fcp_4_taslconnectivity - Fieldtrip-based Connectivity Pipeline, calculate connectivity
% This script will:
%   - load virtual sensor data from fcp_3_beamforming
%   - run time-frequency decomposition methods to extract instantaneous phase at specific
%       frequencies in the signal
%   - estimate functional connectivity using phase-locking metrics (PLV, PLI, wPLI, etc) or
%       amplitude coupling
% Uses essentially the same function from meglegacy
%
% Sonya Bells
% Based on a sctip written by Simeon Wong
% 
% 2016 May 5

%%
config = load_config(paths, paths.name);
config = config.config;
step = 'fcp4';
subj_ds = load_participants(paths,step);
pids = readtable(paths.all_subj_pids);
[subj_match, failure] = ds_pid_match(paths,step);
ssSubjPath = @(x) paths.(subj_match.pid{x});

if strcmp(config.connectivity.method,'pli')
    connfn     = @(H1,H2) abs(mean(sign(angle(H1)-angle(H2)), 1));
elseif strcmp(config.connectivity.method,'plv')
    connfn = @(H1,H2) abs(mean(exp(1i .* (angle(H1)-angle(H2))), 1));
end

%z-score
prefilter = @(ts) zscore(ts);
% filter function
filterfn = @(b,a,ts) filtfilt(b,a,ts);

%% Check data
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

%% Make filters
%load(p.paths.vs(p.subj{1}, p.condition{1}), 'srate');
load([ssSubjPath(1) '/AAL_beamforming_results'], '-mat') %SB - add cond

maxn = 0;
for fq = 1:length(config.connectivity.filt_freqs)
  % calculate filter coefficients
  filtband = config.connectivity.filt_freqs(fq,:);
  
  % Estimate filter order
  transition = mean(filtband) * 0.2;
  freqs = [filtband(1) - transition, filtband(1), filtband(2), filtband(2) + transition];
  fir_n = kaiserord(freqs, [0 1 0], [0.1 0.05 0.1], srate);
  
  maxn = max(maxn, fir_n);
  
  % Make filter coefficients
  fir_coef{fq} = fir1(fir_n, filtband*2/srate, 'bandpass');

%   if (max(length(fir_coef{fq})-1,length(1)-1)*3) > length(ts)
%       fprintf('Warning: For fiter [%d %d] - zero phase filtering can not be used.\n' ,p.filt_freqs(fq,:)) 
%   end
  
end

fprintf('Max filter length: %d samples = %.4f sec.\n', maxn, maxn/srate);



%% Run connectivity
% calculate connectivity over time throughout the x second epoch
band_names = ["theta", "alpha", "beta", "lowgamma", "highgamma"];
all_adjmat = nan(90, 90, length(subj_match.ds), length(config.connectivity.filt_freqs));

%for cc = 1:length(p.condition)
  % fprintf('Running condition %s... \n', p.condition{cc});
  for ss = 1:length(subj_match.ds)
    fprintf('\n\n==================================\nSUBJECT: %s\n', subj_match.pid{ss});

    fprintf('Loading participant %s... \n',subj_match.pid{ss});
    % load virtual sensor data
    %load(p.paths.vs(p.subj{ss}, p.condition{cc}), 'catmatrix', 'srate');
    load([ssSubjPath(ss) '/AAL_beamforming_results'], '-mat'); %SB - add cond
    
    % define some dimensions
    num_samples = size(catmatrix, 1);
    num_trials  = size(catmatrix, 2);
    num_sources = size(catmatrix, 3);
    
    % initialize the individual-subject adjacency matrices
    %  - dimensions [sources] x [sources] x [trials] x [frequencies]
    %adjmat = nan(num_sources, num_sources, num_samples,length(p.filt_freqs)); %orig
    adjmat = nan(num_sources, num_sources, num_trials,length(config.connectivity.filt_freqs)); % in resting state
    
    % loop over frequencies
    %parfor fq = 1:length(p.filt_freqs)
    for fq = 1:length(config.connectivity.filt_freqs) 
      % initialize the temporary adjmat
      fprintf('Analyzing the %s band!\n', band_names(fq)); 
      p_adjmat = nan(num_sources, num_sources, num_trials); % same as RestingSate     
      H_data = complex(zeros(num_samples, num_trials, num_sources)); 
      for tt = 1:num_trials
        fprintf('Working on trial %d... \n',tt);
        % filter all the sources
        for kk = 1:num_sources
          % mean center or z-score the timeseries before filtering
          ts = prefilter(catmatrix(:,tt,kk));          
          % filter the data, calculate the hilbert transform, get the instantaneous phase
          if (max(length(fir_coef{fq})-1,length(1)-1)*3) < length(ts)          
            H_data(:,tt,kk) = hilbert(filtfilt(fir_coef{fq}, 1, ts));
          else
            H_data(:,tt,kk) = hilbert(filter(fir_coef{fq}, 1, ts));
          end
        end
      end
      % calculate connectivity
      for aa = 1:num_sources
        for bb = aa+1:num_sources
          p_adjmat(aa,bb,:) = connfn(H_data(:,:,aa), H_data(:,:,bb));
          p_adjmat(bb,aa,:) = p_adjmat(aa,bb,:);
        end
      end      
      % copy temporary matrix data into subject-matrix
      adjmat(:,:,:,fq) = p_adjmat;
    end

    %   - adjmat is [sources] x [sources] x [time] x [frequencies]
    fprintf('Saving the mat file...');
    save([ssSubjPath(ss) '/fcp_5_adjmat.mat'],'adjmat','-mat','-v7.3')
    all_adjmat(:,:,ss,:) = squeeze(mean(adjmat, 3));
  end
save([paths.anout_grp '/fcp_5_allParticipants_adjmats.mat'],'all_adjmat','-v7.3');

% end
 fprintf('\n\n============== Finished Processing ====================\n');  



