function [ppt_marker_summary, cluster] = getMarkerSummary(paths)

config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp1';

% check for matched MRI and MEG data
subj_match = ds_pid_match(paths,step);


%% look across all participants to get master list of marker types

parfor ss = 1:length(subj_match.ds) % for each participant

    cfg             = [];
    cfg.dataset     = [paths.rawdata '/' subj_match.ds{ss}]; 
    cfg.trialfun    = config.taskFunc; 
    cfg.trialdef    = config.task.trialdef;
    cfg.continuous  = 'yes';

    % formatting in case the field entry is a struct
    switch class(cfg.trialdef.details)
        case 'struct'
            cfg.trialdef.details = struct2table(cfg.trialdef.details, 'AsArray', true);
    end
    
    % extract list of events from the MEG dataset
    eventslist  = ft_read_event(cfg.dataset);

    % get list of possible event types
    uniquetypes{ss} = sort(string(unique({eventslist.type})))';
end

% get unique list of markers across all participants
all_types = [];
for ss = 1:length(subj_match.ds)
    all_types = unique([all_types; uniquetypes{ss}]);
end
% Conclusion: there are 22 unique markers

%% break down by participant
% make readable the unique types for each participant
num_markers = length(all_types);
for m = 1:num_markers
    does_it_exist = cellfun(@(x) any(strcmp(all_types(m), x)), uniquetypes, 'UniformOutput', false);
    ppt_marker_summary(m,:) = cell2mat(does_it_exist);
end

ppt_marker_summary = array2table(ppt_marker_summary);
ppt_marker_summary.Properties.RowNames = all_types;
ppt_marker_summary.Properties.VariableNames = subj_match.pid;
% saved to megneto/analysis

%% which markers are redundant?
    
parfor ss = 1:length(subj_match.ds) % for each participant

    cfg             = [];
    cfg.dataset     = [paths.rawdata '/' subj_match.ds{ss}]; 
    cfg.trialfun    = config.taskFunc; 
    cfg.trialdef    = config.task.trialdef;
    cfg.continuous  = 'yes';

    % formatting in case the field entry is a struct
    switch class(cfg.trialdef.details)
        case 'struct'
            cfg.trialdef.details = struct2table(cfg.trialdef.details, 'AsArray', true);
    end
    
    % extract list of events from the MEG dataset
    eventslist  = ft_read_event(cfg.dataset);

    these_events = string({eventslist.type});
    these_samples = cell2mat({eventslist.sample});
    unique_samples = unique(these_samples);
    
    for samp = 1:length(unique_samples)
        these_markers = these_events(these_samples == unique_samples(samp));
        cluster{ss}(1:length(these_markers), samp) = these_markers;
    end
end
% this produces a cell array called cluster
% there's an entry for each participant
% within each cell for a participant, there is a cluster x sample array
    
    
    
    
    
    
    