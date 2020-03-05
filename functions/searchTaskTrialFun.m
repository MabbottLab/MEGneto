function [trl, eventslist] = searchTaskTrialFun(cfg)

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
%   cfg.trialdef.details.exclude    - what other trigger to exclude
%   cfg.traildef.tEpoch     - time window  eg. [-1.5 1.5]
%
% OUTPUTS:
%   trl
%   eventslist
%
% See also: FCP_1_TASKEPOCHING, FT_DEFINETRIAL
%
% Last updated by: Julie Tseng, 2020-01-07
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.

%% CHECK INPUTS

% formatting in case the field entry is a struct
switch class(cfg.trialdef.details)
    case 'struct'
        cfg.trialdef.details = struct2table(cfg.trialdef.details, 'AsArray', true);
end

% check the number of conditions that are being run
for tt = 1:height(cfg.trialdef.details) % for each condition
    if isempty(cfg.trialdef.details.include{tt}) % if no marker designated
        cfg.trialdef.details.include{tt} = {};
    end
    % I'm not really sure what the heck this is accomplishing? 
    % this is length of the include field, which is a char array aka not a
    % measure of how many there are?
    if length(cfg.trialdef.details.include{tt}) < 1 % if multiple markers are designated?
        warning('This trialdef function does not support multiple included markers per trial. Use another function or write your own.')
    end
    if isempty(cfg.trialdef.details.includeOnce{tt}) 
        cfg.trialdef.details.includeOnce{tt} = {};
    end
    if isempty(cfg.trialdef.details.exclude{tt})
        cfg.trialdef.details.exclude{tt} = {};
    end
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
t0markernum   = find(strcmpi(uniquetypes, cfg.trialdef.markers.t0marker),1);

% setup unfiiltered trials table for all trials
unfiltered_trials               = table(); % instantiate table
unfiltered_trials.t0sample      = [eventslist(eventslistnum == t0markernum).sample]'; % get the sample of each t0marker

% prepare event and event timing columns
unfiltered_trials.event         = cell(height(unfiltered_trials),1); 
unfiltered_trials.eventTiming   = cell(height(unfiltered_trials),1);

%% ACTUAL PROCESSING OF EVENTS

%%% IDENTIFY START AND END SEARCH INTERVAL --------------------------------
for rr = 1:height(unfiltered_trials)
    
    % get the time interval for this trial
    sStart      = unfiltered_trials.t0sample(rr)-1; % from one sample prior to the t0marker
    if rr < height(unfiltered_trials) % if we are not at the last trial
        sEnd    = unfiltered_trials.t0sample(rr+1); % set sEnd to be the starting sample of the next trial
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
marker_present = cellfun(@(x) contains(x,cfg.trialdef.markers.Correct), ...             % return logical vector per trial of time interval
                        unfiltered_trials.event, 'UniformOutput', false);
marker_present = cellfun(@(x) find(x,1), marker_present, 'UniformOutput', false);       % find the 1 if it exists and determine the sample
marker_present(cellfun(@isempty,marker_present)) = {0};                                 % handle empties
marker_present = cell2mat(marker_present);       % make it into a matrix
deleteRows = marker_present == 0;
marker_present = marker_present(marker_present > 0);  

%% ASSEMBLE TRL MATRIX

%%% INITIALIZE VARIABLES --------------------------------------------------
trl         = []; % the trial matrix
trl_counts  = zeros(height(cfg.trialdef.details));   % track num trials

%%% FILTER THE TRIAL LIST FOR ONES THAT DO CONTAIN DESIGNATED MARKER-------
selected_trials                 = unfiltered_trials;
selected_trials(deleteRows,:)   = [];

% remove or modify once you set up multi-condition
tt = 1;

%%% FORM TRL MATRIX -------------------------------------------------------
if ~cfg.trialdef.details.countOnly % if you don't want only counts
    
    % trl_tmp: 
    %       - COL 1 and 2:  samples of [-2, 2] around t0marker + t0shift
    %       - COL 3:        time in msec of lower bound from t = 0
    %       - COL 4:        don't worry, tracking num_conditions but here = 1
    trl_tmp         = zeros(height(selected_trials), 4 + length(cfg.trialdef.details.include));
    trl_tmp(:,1:2)  = selected_trials.t0sample ... % around the sample of t0marker
                        + timeToSamp(cfg.trialdef.parameters.tEpoch) ... % get [-2, 2]
                        + timeToSamp(cfg.trialdef.parameters.t0shift); % shift by t0shift
    trl_tmp(:,3)    = timeToSamp(cfg.trialdef.parameters.tEpoch(1));
    trl_tmp(:,4)    = tt;

    % combine include once and includes for stats
    includes        = reshape(cat(2, cfg.trialdef.details.include(tt), cfg.trialdef.details.includeOnce(tt)), [], 1);
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
    
    % throw out if they are too close to the beginning or end
    trl_tmp = trl_tmp(trl_tmp(:,1) >= (1+timeToSamp(0.6)) & trl_tmp(:,2) <= hdr.nSamples, :);
    trl = [trl; trl_tmp];
end
trl_counts(tt) = height(selected_trials); 
    
%% MAKE SOME STATEMENTS ABOUT PROGRESS/CONCLUSIONS
fprintf('======================================\n');
fprintf('============ TRIAL COUNTS ============\n');
fprintf('======================================\n');

for tt = 1:height(cfg.trialdef.details)
    fprintf('%d\t- %s\n', trl_counts(tt), cfg.trialdef.details.name{tt});
end
fprintf('--------------------------------------\n');
fprintf('%d\t- TOTAL\n', height(unfiltered_trials));
fprintf('--------------------------------------\n');

fprintf('\n\n');

% sort trl matrix
trl = sortrows(trl);

end
