function [trl,eventslist]= taskTrialFun(cfg)
% Defines trials with t = 0 at a stimulus onset marker, based on  motor
% response
% Sonya Bells
% Based on Simeon Wong's code (includesBasedTrialFun)
%%Required Inputs %%%
% cfg.dataset             - path to CTF dataset
% cfg.trialdef.t0marker   - defines the t = 0 of each trial (i.e. presentation trigger)
% cfg.trialdef.parameters.t0shift    - time in seconds to offset t0marker (presentation delay)
% cfg.trialdef.details.trigger    - name of trigger
% cfg.trialdef.details.include    - what other trigger to include
% cfg.trialdef.details.exclude    - what other trigger to exinclude
%cfg.trialdef.parameters.tSearch    - time interval in seconds to search for inclusion/exclusion markers as
%                           1x2 vector   eg. [-0.5, 1.0]
% cfg.traildef.tEpoch     - time window  eg. [-1.5 1.5]
%%
% check inputs
cfg.trialdef.details = struct2table(cfg.trialdef.details, 'AsArray', true);
for tt = 1:height(cfg.trialdef.details)
    if isempty(cfg.trialdef.details.include{tt})
        cfg.trialdef.details.include{tt} = {};
    end
    if isempty(cfg.trialdef.details.includeOnce{tt})
        cfg.trialdef.details.includeOnce{tt} = {};
    end
    if isempty(cfg.trialdef.details.exclude{tt})
        cfg.trialdef.details.exclude{tt} = {};
    end
end

if ~isfield(cfg.trialdef.markers, 't0marker')
    error('Missing t0marker index');
end

if ~isfield(cfg.trialdef.parameters, 'tSearch') %|| any(size(cfg.trialdef.parameters.tSearch) ~= [1 2])
    error('Missing or invalid tSearch');
end

if ~isfield(cfg.trialdef.parameters, 'tEpoch') %|| any(size(cfg.trialdef.parameters.tEpoch) ~= [1 2])
    error('Missing or invalid tEpoch');
end


%% load metadata from dataset
%%% meg dataset header
% mainly just for sampling rate
hdr = ft_read_header(cfg.dataset);

timeToSamp = @(t) round(t * hdr.Fs);
sampToTime = @(s) s / hdr.Fs;

%%% load events
% extract list of events from the MEG dataset
eventslist= ft_read_event(cfg.dataset);

% get list of possible event types
uniquetypes = unique({eventslist.type});

% assign index to events (instead of just string)
eventslistnum = cellfun(@(str) find(strcmp(uniquetypes, str),1), {eventslist.type});


%% process events

% identify all the t=0 stimulus markers, and search for other markers within the search interval
t0markernum = find(strcmpi(uniquetypes, cfg.trialdef.markers.t0marker),1);


unfiltered_trials = table();
unfiltered_trials.t0sample = [eventslist(eventslistnum == t0markernum).sample]';
unfiltered_trials.event = cell(height(unfiltered_trials),1);
unfiltered_trials.eventTiming = cell(height(unfiltered_trials),1);

for rr = 1:height(unfiltered_trials)
    sStart = unfiltered_trials.t0sample(rr) + timeToSamp(cfg.trialdef.parameters.tSearch(1));
    sEnd = unfiltered_trials.t0sample(rr) + timeToSamp(cfg.trialdef.parameters.tSearch(2));
    unfiltered_trials.event{rr} = {eventslist([eventslist.sample] >= sStart & [eventslist.sample] <= sEnd).type};
    unfiltered_trials.eventTiming{rr} = {eventslist([eventslist.sample] >= sStart & [eventslist.sample] <= sEnd).sample};
end


% loop through each of the trial requirements given in trialdef.details
% create a list of trials that fit that requirement
filtered_trial_list = false(height(unfiltered_trials), height(cfg.trialdef.details));

for tt = 1:height(cfg.trialdef.details)
    for rr = 1:height(unfiltered_trials)
        if ~isempty(cfg.trialdef.details.include{tt})
            trialIncludes = cellfun(@(x) any(strcmpi(unfiltered_trials.event{rr}, x)), cfg.trialdef.details.include(tt));
        else
            trialIncludes = true;
        end
        if ~isempty(cfg.trialdef.details.includeOnce{tt})
            trialIncludeOnce = cellfun(@(x) sum(strcmpi(unfiltered_trials.event{rr}, x)) == 1, cfg.trialdef.details.includeOnce(tt));
        else
            trialIncludeOnce = true;
        end
        if ~isempty(cfg.trialdef.details.exclude{tt})
            trialExcludes = cellfun(@(x) any(strcmpi(unfiltered_trials.event{rr}, x)), cfg.trialdef.details.exclude(tt));
        else
            trialExcludes = false;
        end
        %disp(trialIncludeOnce)
        filtered_trial_list(rr,tt) = all(trialIncludes) & all(trialIncludeOnce) & ~any(trialExcludes); %all(trialIncludes) & all(trialIncludeOnce) & ~any(trialExcludes);
    end
end

% disp(filtered_trial_list(:,1)');
% double check that there aren't any trials that are double-classified
cfg.trialdef.details.countOnly = logical(cfg.trialdef.details.countOnly);
filtered_trial_keeponly = filtered_trial_list(:,~cfg.trialdef.details.countOnly);
filtered_trial_totalclasses = sum(filtered_trial_keeponly,2);
idx = find(filtered_trial_totalclasses >= 2);

for kk = 1:length(idx)
    warning('Trial %d belongs to more than one classification!\n', idx(kk));
end

%% assemble trl matrix
trl = [];
trl_counts = zeros(height(cfg.trialdef.details));
% disp((cfg.trialdef.details));


for tt = 1:height(cfg.trialdef.details)
    
    if height(cfg.trialdef.details) == 1
        
        deleteRows = find(filtered_trial_list(:,tt) == 0)';
        selected_trials = unfiltered_trials; %(filtered_trial_list(:,tt),:,:);
        selected_trials(deleteRows,:) = [];
    else
        selected_trials = unfiltered_trials(filtered_trial_list(:,tt),:,:);
    end
    
    if ~cfg.trialdef.details.countOnly(tt)
        trl_tmp = zeros(height(selected_trials), 4 + length(cfg.trialdef.details.include{tt}));
        trl_tmp(:,1) = selected_trials.t0sample + timeToSamp(cfg.trialdef.parameters.tEpoch(1));
        trl_tmp(:,2) = selected_trials.t0sample + timeToSamp(cfg.trialdef.parameters.tEpoch(2));
        trl_tmp(:,3) = timeToSamp(cfg.trialdef.parameters.tEpoch(1));
        trl_tmp(:,4) = tt;
        
        % shift all trials by t0shift
        trl_tmp(:,1:2) = trl_tmp(:,1:2) + timeToSamp(cfg.trialdef.parameters.t0shift);
        
        % combine includeonce and includes for stats
        includes = reshape(cat(2, cfg.trialdef.details.include(tt), cfg.trialdef.details.includeOnce(tt)), [], 1);
        includes = includes(~isempty(includes));
        % loop through trials and pull latency
        for rr = 1:size(trl_tmp,1)
            % get latencies
            for mm = 1:length(includes)
                disp(mm)
                trl_tmp(rr,4+mm) = selected_trials.eventTiming{rr}{strcmpi(selected_trials.event{rr}, includes{mm})};
            end
            
            % convert to milliseconds wrt t=0
            trl_tmp(rr,5:end) = sampToTime(trl_tmp(rr,5:end) - selected_trials.t0sample(rr));
        end
        
        % drop trials that are outside the end of the file
        trl_tmp = trl_tmp(trl_tmp(:,1) >= 1 & trl_tmp(:,2) <= hdr.nSamples, :);
        
        trl = [trl; trl_tmp];
    end
    trl_counts(tt) = height(selected_trials);
    
end



fprintf('======================================\n');
fprintf('============ TRIAL COUNTS ============\n');
fprintf('======================================\n');

for tt = 1:height(cfg.trialdef.details)
    fprintf('%d\t- %s\n', trl_counts(tt), cfg.trialdef.details.name{tt});
end
fprintf('--------------------------------------\n');
fprintf('%d\t- TOTAL\n', height(unfiltered_trials));
fprintf('--------------------------------------\n');
fprintf('Average ISI is %.4f sec', sampToTime(mean(diff(unfiltered_trials.t0sample(any(filtered_trial_keeponly,2))))));

fprintf('\n\n');


% sort trl matrix
trl = sortrows(trl);

end




