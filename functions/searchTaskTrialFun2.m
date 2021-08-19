function [trl,eventslist]= searchTaskTrialFun2(cfg)
% Defines trials with t = 0 at a stimulus onset marker, based on a motor
% response. Defines trial search window as from one stim onset marker to
% the next or EOF otherwise.
% Ming Scott's code
% Adapted from
% Sonya Bells's code (taskTrialFun)
% Based on Simeon Wong's code (includesBasedTrialFun)
%%Required Inputs %%%
% cfg.dataset             - path to CTF dataset
% cfg.trialdef.t0marker   - defines the t = 0 of each trial (i.e. presentation trigger)
% cfg.trialdef.parameters.t0shift    - time in seconds to offset t0marker (presentation delay)
% cfg.trialdef.details.trigger    - name of trigger
% cfg.trialdef.details.include    - what other trigger to include
% cfg.trialdef.details.exclude    - what other trigger to exinclude
% cfg.traildef.tEpoch     - time window  eg. [-1.5 1.5]
%%
% check inputs
switch class(cfg.trialdef.details)
    case 'struct'
        cfg.trialdef.details = struct2table(cfg.trialdef.details, 'AsArray', true);
end
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
    sStart = unfiltered_trials.t0sample(rr)-1;
    if rr < height(unfiltered_trials)
        sEnd = unfiltered_trials.t0sample(rr+1);
    else
        sEnd = eventslist(end).sample;
    end %Defines marker search interval as from one sample before t0marker
% to next t0marker
    unfiltered_trials.event{rr} = {eventslist([eventslist.sample] >= sStart & [eventslist.sample] <= sEnd).type};
    unfiltered_trials.eventTiming{rr} = {eventslist([eventslist.sample] >= sStart & [eventslist.sample] <= sEnd).sample};
end
% for ii = 1:height(unfiltered_trials)
%     if unfiltered_trials.t0sample(ii) + timeToSamp(cfg.trialdef.parameters.t0shift) + timeToSamp(cfg.trialdef.parameters.tEpoch(1)) - timeToSamp(0.1) < 0 %this is hard coded!!
%         unfiltered_trials(ii,:) = [];
%     end
% end
% This line checks if any of the markers designated for trial inclusion are
% present within each trial.
contains = cellfun(@(x) cellfun(@(y) ...
        strcmp(x,y),...
        unfiltered_trials.event,'UniformOutput',false),... x
        cfg.trialdef.details.include,'UniformOutput',false); % y
lenc = length(contains);
all = cell(length(contains),1);
for i=1:length(contains{1})
    tmp = zeros(length(contains{1}{i}), 1)';
	for ii=1:lenc
        tmp = [tmp; contains{ii}{i}];
    end
    all{i} = sum(tmp);
end
contains = all;
% Finds the first instance of a marker designated for inclusion in each
% trial, else returns and empty string
contains = cellfun(@(x) find(x,1),contains,'UniformOutput',false);
% Replaces the empty strings with zeros
contains(cellfun(@isempty,contains)) = {0};
% Converts to a numerical vector, where nonzero integer values are truthy
filtered_trial_list = cell2mat(contains);

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
            for mm = 1:length(selected_trials.event{rr})
                tmp = strfind(selected_trials.event{rr}{mm}, 'Correct');
                if ~isempty(tmp)
                    trl_tmp(rr,5) = selected_trials.eventTiming{rr}{strfind(selected_trials.event{rr}{mm}, 'Correct')};
                end
            end
            
            % convert to milliseconds wrt t=0
            trl_tmp(rr,5:end) = sampToTime(trl_tmp(rr,5:end) - selected_trials.t0sample(rr));
        end
        
        % drop trials that are outside the end of the file
        trl_tmp = trl_tmp(trl_tmp(:,1) >= (1+timeToSamp(0.1)+timeToSamp(0.5)) & trl_tmp(:,2) <= hdr.nSamples, :);
        
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
% fprintf('Average ISI is %.4f sec', sampToTime(mean(diff(unfiltered_trials.t0sample(any(filtered_trial_keeponly,2))))));

fprintf('\n\n');


% sort trl matrix
trl = sortrows(trl);

end
