function [trl, eventslist] = trialfun_MSIT_correct(cfg)

% SEARCHTASKTRIALFUN is the function that filters and organizes trials by 
% the designated markers. This function is called by:
%   main > fcp_1_taskepoching > ft_definetrial
%
% INPUTS:
%   cfg.dataset             - path to CTF dataset
%   cfg.trialdef.t0marker   - defines the t = 0 of each trial (i.e. presentation trigger)
%   cfg.trialdef.parameters.t0shift    - time in seconds to offset t0marker (presentation delay)
%   cfg.trialdef.details.trigger    - name of trigger
%   cfg.trialdef.details.include    - what other trigger to include
%   cfg.traildef.tEpoch     - time window  eg. [-1.5 1.5]
%
% OUTPUTS:
%   trl
%   eventslist
%

%% LOAD METADATA

% here's the header
hdr         = ft_read_header(cfg.dataset);

% create functions to calculate:
timeToSamp  = @(t) round(t * hdr.Fs); % time x sampling rate = sample
sampToTime  = @(s) s / hdr.Fs; % sample / sampling rate = time

% extract list of events from the MEG dataset
eventslist  = ft_read_event(cfg.dataset);

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
deleteRows = marker_present == 0; % get rows of incorrect trials
marker_present = marker_present(marker_present > 0);  

%% ASSEMBLE TRL MATRIX

%%% INITIALIZE VARIABLES --------------------------------------------------
trl         = []; % the trial matrix

%%% FILTER THE TRIAL LIST FOR ONES THAT DO CONTAIN DESIGNATED MARKER-------
selected_trials                 = unfiltered_trials;
selected_trials(deleteRows,:)   = [];

%%% FORM TRL MATRIX -------------------------------------------------------
for rr = 1:size(selected_trials, 1)
    t0_index = contains(selected_trials{rr, 'event'}{1}, 'Button');
    these_eventTimings = cell2mat(selected_trials{rr, 'eventTiming'}{1});
    this_t0 = these_eventTimings(t0_index);
    
    trl(rr,1) = this_t0 - timeToSamp(cfg.trialdef.prestim);
    trl(rr,2) = this_t0 + timeToSamp(cfg.trialdef.poststim);
    trl(rr,3) = -1 * timeToSamp(cfg.trialdef.prestim);
    trl(rr,4) = string(selected_trials{rr, 'condition'}{1}) == "CONG";
end
    
%% print some messages

disp(sprintf("Started with %d trials,", size(unfiltered_trials,1)));
disp(sprintf("ended up with %d correct response trials.", size(selected_trials,1)));
disp("Fourth column of trl matrix is 1 if CONGRUENT, 0 if INCONGRUENT.")

end
