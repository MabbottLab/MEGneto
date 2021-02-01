function [reactionTimeTable] = reactionTimes(cfg)

% reactionTimes is a function that extracts participant reaction
% times. That is, the time between stimulus presentation and button press.
% Currently this function is used for the Speed of Thinking TP & EL 
% dataset.

% INPUTS:
%   cfg.dataset                     - path to CTF dataset
%   cfg.trialdef.t0marker           - defines the t = 0 of each trial 
%                                     (i.e. presentation trigger)
%   cfg.trialdef.parameters.t0shift - time in seconds to offset t0marker 
%                                     (presentation delay)
%   cfg.trialdef.details.trigger    - name of trigger
%   cfg.trialdef.details.include    - what other trigger to include
%   cfg.trialdef.details.exclude    - what other trigger to exclude
%   cfg.traildef.tEpoch             - time window  eg. [-1.5 1.5]
%
% OUTPUTS: reactionTimeTable 
%               reactionTimeTable.t0sample       - sample number of stimulus 
%                                                  presentation
%               reactionTimeTable.event          - marker label for the  
%                                                  event of interest
%               reactionTimeTable.eventTiming    - sample number of event  
%                                                  of interest
%               reactionTime.Table.reactionTimes - time (seconds) between
%                                                  stimulus presentation 
%                                                  and button press
%
% Note: to convert sample number to time, simply divide the sample number
% by the sampling frequency. Similarly, to convert time to sample number,
% multiple the time by the sampling frequency. 
%
% Last updated by: Dunja Matic, 2021-01-29

%% CHECK INPUTS

% formatting in case the field entry is a struct
switch class(cfg.trialdef.details)
    case 'struct'
        cfg.trialdef.details = struct2table(cfg.trialdef.details, 'AsArray', true);
end

% check for designated t0 marker label
if ~isfield(cfg.trialdef.markers, 't0marker')
    error('Missing t0marker index');
end

% check that a epoch interval was given
if ~isfield(cfg.trialdef.parameters, 'tEpoch') 
    error('Missing or invalid tEpoch');
end

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
t0markernum   = find(ismember(uniquetypes, cfg.trialdef.markers.t0marker));

% setup unfiiltered trials table for all trials
reactionTimeTable               = table(); % instantiate table
reactionTimeTable.t0sample      = [eventslist(ismember(eventslistnum, t0markernum)).sample]'; % get the sample of each t0marker

% prepare event and event timing columns
reactionTimeTable.event         = cell(height(reactionTimeTable),1); 
reactionTimeTable.eventTiming   = cell(height(reactionTimeTable),1);

%% ACTUAL PROCESSING OF EVENTS

%%% IDENTIFY START AND END SEARCH INTERVAL --------------------------------
for rr = 1:height(reactionTimeTable)
    
    % get the time interval for this trial
    sStart      = reactionTimeTable.t0sample(rr)-1; % from one sample prior to the t0marker
    if rr < height(reactionTimeTable) % if we are not at the last trial
        sEnd    = reactionTimeTable.t0sample(rr+1)-1; % set sEnd to be the sample before the next trial
    else
        sEnd    = eventslist(end).sample; % otherwise, set it to be the last sample
    end 
    
    % get the epoch of events
    reactionTimeTable.event{rr} = {eventslist([eventslist.sample] >= sStart ... for sStart < samples < sEnd
                                & [eventslist.sample] <= sEnd).type}; % get the label of the event
    reactionTimeTable.eventTiming{rr} = {eventslist([eventslist.sample] >= sStart ...
                                      & [eventslist.sample] <= sEnd).sample}; % get the sample ordinal 
end

for i = 1:length(reactionTimeTable.event)
    current_cell = reactionTimeTable.event{i};
    quick_check = contains(current_cell, 'orrect');
    do_any_contain_1 = any(quick_check);
    if do_any_contain_1 == 1
        for j = 1:length(current_cell)
            current_event = current_cell{j};
            if contains(current_event, 'orrect')
                reactionTimeTable.event{i} = current_event;
                reactionTimeTable.eventTiming{i} = reactionTimeTable.eventTiming{i}{j};
            end          
        end
    elseif do_any_contain_1 == 0
        reactionTimeTable.event{i} = NaN;
        reactionTimeTable.eventTiming{i} = NaN;
    end
end

for i = 1:height(reactionTimeTable)
    reactionTimeTable.reactionTimes(i) = sampToTime(reactionTimeTable.eventTiming{i} - reactionTimeTable.t0sample(i));
end 
  
end