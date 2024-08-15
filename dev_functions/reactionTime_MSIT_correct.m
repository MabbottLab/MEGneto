function rt_all = reactionTime_MSIT_correct(cfg)

%% LOAD METADATA

% here's the header
hdr         = ft_read_header(cfg.dataset);

% create functions to calculate:
timeToSamp  = @(t) round(t * hdr.Fs); % time x sampling rate = sample
sampToTime  = @(s) s / hdr.Fs; % sample / sampling rate = time

% extract list of events from the MEG dataset
eventslist  = ft_read_event(cfg.dataset);

% handle FIT_004's missing buttons
if contains(cfg.dataset, "FIT-004-01")
    eventslist = ft_read_event('/media/mabbotthpf/datasets/FitABCS/FIT_004/ses-01/meg/FIT-004-01_RESEARCH_20220811_MSIT.ds/MarkerFile_addButtons_20230420.mat');
end

% get list of possible event types
uniquetypes = unique({eventslist.type});

% assign index to events (instead of just string)
eventslistnum = cellfun(@(str) find(strcmp(uniquetypes, str),1), {eventslist.type});

% identify all the t=0 stimulus markers, and search for other markers within the search interval
tripletnum   = find(contains(uniquetypes, cfg.trialdef.eventtype));

% setup unfiltered trials table for all trials
unfiltered_trials               = table(); % instantiate table
unfiltered_trials.triplet      = [eventslist(ismember(eventslistnum, tripletnum)).sample]'; % get the sample of each t0marker

% prepare event and event timing columns
unfiltered_trials.event         = cell(height(unfiltered_trials),1); 
unfiltered_trials.eventTiming   = cell(height(unfiltered_trials),1);
unfiltered_trials.condition     = {(eventslist(ismember(eventslistnum, tripletnum)).type)}';

%% ACTUAL PROCESSING OF EVENTS

%%% IDENTIFY START AND END SEARCH INTERVAL --------------------------------
for rr = 1:height(unfiltered_trials)
    
    % handle sStart
    sStart = unfiltered_trials.triplet(rr);
    % handle sEnd
    if rr < height(unfiltered_trials) % if we are not at the last trial
        sEnd    = unfiltered_trials.triplet(rr+1)-1; % set sEnd to be the sample before the next trial
    else
        sEnd    = eventslist(end).sample; % otherwise, set it to be the last sample
    end 
    
    % get the epoch of events
    unfiltered_trials.event{rr} = {eventslist([eventslist.sample] >= sStart ... for sStart < samples < sEnd
                                & [eventslist.sample] <= sEnd).type}; % get the label of the event
    unfiltered_trials.eventTiming{rr} = {eventslist([eventslist.sample] >= sStart ...
                                      & [eventslist.sample] <= sEnd).sample}; % get the sample ordinal 
end

%%% DETERMINE DESIGNATED MARKER PRESENCE WITHIN EACH TRIAL ----------------

% marker_present keeps track of the sample at which the designated marker
% occurs for each trial. 
marker_present = cellfun(@(x) any(ismember(x,'Correct')), ...             % return logical vector per trial of time interval
                        unfiltered_trials.event, 'UniformOutput', false);
marker_present = cell2mat(marker_present);       % make it into a matrix

%% ASSEMBLE TRL MATRIX

%%% INITIALIZE VARIABLES --------------------------------------------------
rt_all         = NaN(size(unfiltered_trials,1), 3); % the trial matrix
rt_all(:,3)    = marker_present;

%%% FILTER THE TRIAL LIST FOR ONES THAT DO CONTAIN DESIGNATED MARKER-------
% selected_trials                 = unfiltered_trials;
% selected_trials(deleteRows,:)   = [];

%%% FORM TRL MATRIX -------------------------------------------------------
for rr = 1:size(unfiltered_trials, 1)
    between_indices = [find(contains(unfiltered_trials{rr, 'event'}{1}, 'CONG')) ...
                    find(contains(unfiltered_trials{rr, 'event'}{1}, 'Button'))];
    if length(between_indices) == 2
        these_eventTimings = cell2mat(unfiltered_trials{rr, 'eventTiming'}{1});
        rt = sampToTime(diff(these_eventTimings(between_indices)));
        rt_all(rr,1) = rt(1); % reaction time in seconds
    end 
    rt_all(rr,2) = string(unfiltered_trials{rr, 'condition'}{1}) == "CONG"; % 1 if congruent, 0 if incongruent
end

rt_all = array2table(rt_all, 'VariableNames', {'ReactionTime', 'CONG (1) or INCONG (0)', 'CORRECT (1) or INCORRECT (0)'});

%% print some msgs

fprintf('For ds file: %s\n', cfg.dataset)
fprintf('They got %d out of %d correct (%0.1f percent).\n', sum(rt_all{:,3}), size(rt_all,1), ...
                            100*sum(rt_all{:,3})/size(rt_all,1))
fprintf('Their mean reaction time was %0.1f seconds.\n', nanmean(rt_all{rt_all{:,3} == 1,1}))

end
