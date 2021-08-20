% we have every clip marker as sample number and we want to epoch 400 ms
% after the first clip marker to 50 ms before the second clip marker over
% and over again 

% define time interval for this study
timeEpoch = 0.4 + (sampToTime(t0marker(next))-0.05); 
%%% FORM TRL MATRIX -------------------------------------------------------
% trl_tmp:
%       - COL 1 and 2:  samples of [-2, 2] around t0marker + t0shift
%       - COL 3:        time in msec of lower bound from t = 0
%       - COL 4:        don't worry, tracking num_conditions but here = 1
trl_tmp         = zeros(height(selected_trials), 4 + length(cfg.trialdef.details.include));
next_change = selected_trials.t0sample;

epoch_interval = [cfg.trialdef.parameters.tEpoch(1), nextchange(next)-cfg.trialdef.parameters.tEpoch(2)];
trl_tmp(:,1:2)  = selected_trials.t0sample ... % around the sample of t0marker
    + timeToSamp(cfg.trialdef.parameters.tEpoch) ... % get [-2, 2]
    + timeToSamp(cfg.trialdef.parameters.t0shift); % shift by t0shift
trl_tmp(:,3)    = timeToSamp(cfg.trialdef.parameters.tEpoch(1));
trl_tmp(:,4)    = tt;

% combine include once and includes for stats
includes        = reshape(cat(2, cfg.trialdef.details.include(tt)), [], 1);
includes        = includes(~isempty(includes));

% loop through trials and pull latency
for rr = 1:size(trl_tmp,1) % for each trial
    % get latencies
    for mm = 1:length(includes)
        % trl_tmp(rr,4+mm) = selected_trials.eventTiming{rr}{strcmpi(selected_trials.event{rr}, includes{mm})};
        % place in COL 5 = (sample where designated marker occurred) + (t0shift in samples)
        trl_tmp(rr,4+mm) = selected_trials.eventTiming{rr}{marker_present(rr)} + timeToSamp(cfg.trialdef.parameters.t0shift);
    end
    
    % convert to milliseconds wrt t=0
    trl_tmp(rr,5:end) = sampToTime(trl_tmp(rr,5:end) - selected_trials.t0sample(rr));
end