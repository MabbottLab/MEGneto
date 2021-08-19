function [trl, eventslist] = freeviewingTaskTrialFun(cfg)

% Compared to the standard MEG pipeline, the following lines were edited:
% This function is edited significantly and should be walked through by
% the user line by line if they wish to understand it better in comparison
% to the standard MEG pipeline's taskTrialFun.m.

% FREEVIEWINGTASKTRIALFUN is the function that filters and organizes trials by 
% the designated marker 'clipChange'. This function is called by:
%   main > fcp_1_taskepoching > ft_definetrial
%
% INPUTS:
%   cfg.dataset             - path to CTF dataset
%   cfg.trialfun            - config.taskFunc (specifies to use this freeviewingTaskTrialFun);
%   cfg.trialdef            - config.task.trialdef (specified in settings of JSON config);

% OUTPUTS:
%   trl
%   eventslist
%
% See also: FCP_1_TASKEPOCHING, FT_DEFINETRIAL
%
% Last updated by: Dunja Matic, 2021-08-13
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

% extract current participant (since we extract markers one at ppt at a
% time)
pos = find(cfg.dataset == '/', 1, 'last');
ppt = extractAfter(cfg.dataset, pos);

% load clip markers file - will be stored in clipMarkers_allPpts variable
clipMarkers_allPpts = load([clipTimes_path '/clipMarkers_allPpts.mat']);
clipMarkers_allPpts = clipMarkers_allPpts.clipMarkers_allPpts;

% create eventslist structure
current_row = clipMarkers_allPpts(strcmp(clipMarkers_allPpts.Properties.RowNames, ppt), :);
current_col = current_row{:,3};
eventslist = struct();
movieOrder = current_row{:,2}{1};
n = 1;

% assign a movie number to each sample to facilitate rearrangment later on
for i = 1:length(current_col{1})
    for j = 1:length(current_col{1}{i})
        eventslist(n).sample = current_col{1}{i}(j);
        eventslist(n).type = 'clipChange';
        eventslist(n).movie = movieOrder(i); % add movie number for each sample
        n = n + 1;
    end
end

% get list of possible event types
uniquetypes = unique({eventslist.type});

% assign index to events (instead of just string)
eventslistnum = cellfun(@(str) find(strcmp(uniquetypes, str),1), {eventslist.type});

% identify all the t=0 stimulus markers, and search for other markers within the search interval
t0markernum   = find(ismember(uniquetypes, cfg.trialdef.markers.t0marker));

% setup unfiiltered trials table for all trials
unfiltered_trials               = table(); % instantiate table
unfiltered_trials.t0sample      = [eventslist(ismember(eventslistnum, t0markernum)).sample]'; % get the sample of each t0marker
unfiltered_trials.movieNumber   = [eventslist(ismember(eventslistnum, t0markernum)).movie]'; % get the movie of each t0marker

% prepare event and event timing columns
unfiltered_trials.event         = cell(height(unfiltered_trials),1); 
unfiltered_trials.eventTiming   = cell(height(unfiltered_trials),1);

%% ACTUAL PROCESSING OF EVENTS

%%% IDENTIFY START AND END SEARCH INTERVAL --------------------------------
for rr = 1:height(unfiltered_trials)
    
    % get the time interval for this trial
    sStart      = unfiltered_trials.t0sample(rr)-1; % from one sample prior to the t0marker
    if rr < height(unfiltered_trials) % if we are not at the last trial
        sEnd    = unfiltered_trials.t0sample(rr+1)-1; % set sEnd to be the sample before the next trial
    else
        sEnd    = eventslist(end).sample; % otherwise, set it to be the last sample
    end 
    
    % get the epoch of events
    unfiltered_trials.event{rr} = {eventslist([eventslist.sample] >= sStart ... for sStart < samples < sEnd
                                & [eventslist.sample] <= sEnd).type}; % get the label of the event
    unfiltered_trials.eventTiming{rr} = {eventslist([eventslist.sample] >= sStart ...
                                      & [eventslist.sample] <= sEnd).sample}; % get the sample ordinal 
end

% sort the rows in ascending order such that all runs are ordered into 
% movies 1-5 or 6-10 (this way we can collapse across trials and know clip
% content is consistent in each trial)
unfiltered_trials = sortrows(unfiltered_trials, [2 1]);

%%% DETERMINE DESIGNATED MARKER PRESENCE WITHIN EACH TRIAL ----------------

% marker_present keeps track of the sample at which the designated marker
% occurs for each trial. 
marker_present = cellfun(@(x) any(ismember(x,cfg.trialdef.markers.Correct)), ...             % return logical vector per trial of time interval
                        unfiltered_trials.event, 'UniformOutput', false);
marker_present = cell2mat(marker_present);       % make it into a matrix
deleteRows = marker_present == 0; % get rows of incorrect trials
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

% throw out if they are too close to the beginning or end
%trl_tmp = trl_tmp(trl_tmp(:,1) >= (1+timeToSamp(0.6)) & trl_tmp(:,2) <= hdr.nSamples, :);
trl = [trl; trl_tmp];

trl_counts(tt) = height(selected_trials); 
end
