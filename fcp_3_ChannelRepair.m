function fcp_3_ChannelRepair(paths)

% FCP_3_CHANNELREPAIR remove and repair bad channels detected from fcp_1, 
% but that we held off on removing until the data had been ICA-cleaned. 
% 
% NOTES:
%   - Ensure that subj_fcp3.csv is populated with the subject IDs of
%   participants you want to include after checking over initial results. 
%
% INPUTS:
%   paths               =   struct defining paths to data, participant
%                           folders, analysis folders, config files, etc. 
%
% OUTPUTS:
%   ft_meg_fullyProcessed.mat
%       Head motion removed, muscle + jump artifacts removed, 3rd order
%       gradients accounted for, ICA-cleaned, bad channels repaired data. 
%
% See also: 
%
% Last updated by: Julie Tseng, 2020-01-08
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.

%% SET UP LOGGING FILE

right_now = clock;
log_filename = [paths.conf_dir '/log_' sprintf('%d%d%d', right_now(1:3))];
diary(log_filename)

fprintf('\n\n%d:%d:%02.f       Now running **%s**.\n', ...
    right_now(4:6), mfilename)

%% SETUP

% load config JSON with analysis parameters
config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp3';

% check for matched MRI and MEG data
subj_match  = ds_pid_match(paths,step);
ssSubjPath  = @(x) paths.(subj_match.pid{x});
if isempty(subj_match) % if there are no full sets of data
    error('No participants selected')
end

% save the final participants list to the subj_match_fcp1 CSV
write_match_if_not_empty(paths,step);

% check for multiple *.ds folders for participants
if length(unique(subj_match.pid)) ~= length(subj_match.pid)
    error('More than one ds per participant!')
end

% pull in fcp_2/2_5 output, containing ICA-cleaned data
fcp2_output = loadjson([paths.anout_grp '/fcp2_5_output.json']);
fcp2_output = recursive_json_struct_string_to_func(fcp2_output);

%% REPAIR BAD CHANNELS

rangeOFsubj = 1:height(subj_match);

disp('Starting channel repair...');
for ss = rangeOFsubj
    right_now = clock;
    fprintf('%d:%d:%02.f       Working on SUBJECT: %s!\n', ...
        right_now(4:6), subj_match.pid{ss})

%%% LOAD DATA -------------------------------------------------------------
    load([ssSubjPath(ss) '/' fcp2_output.preprocessedData_cfg],'-mat','data');

%%% IF THERE ARE BAD CHANNELS TO REMOVE
    % get info about bad channels
    channel_check = dir([paths.(subj_match.pid{ss}) '/badChannels.json']);
    if config.cleaningOptions.rmBadChannels == 1 && channel_check.bytes > 5

    %%% FIND NEIGHBOURS ---------------------------------------------------
        disp('Finding neighbours ...')

        cfg               = [];
        cfg.method        = 'distance'; %'distance', 'triangulation' or 'template'
        cfg.neighbourdist = 5; %number, maximum distance between neighbouring sensors (only for 'distance')
        cfg.template      = 'ctf151_neighb.mat'; % name of the template file, e.g. CTF275_neighb.mat
        neighbours        = ft_prepare_neighbours(cfg, data);

    %%% REPAIR BAD CHANNELS -----------------------------------------------
        disp('Repairing bad channels ...')

        % run repair
        cfg             = [];
        cfg.method      = 'weighted';        %'average', 'spline' or 'slap' (default = 'weighted')
        cfg.badchannel  = loadjson([paths.(subj_match.pid{ss}) '/badChannels.json']); 
        cfg.neighbours  = neighbours; %neighbourhood structure, see also FT_PREPARE_NEIGHBOURS
        cfg.senstype    = 'meg';
        data            = ft_channelrepair(cfg, data);

        % save result
        save([ssSubjPath(ss) '/ft_meg_fullyProcessed'],'data','-v7.3')
        disp('Done Removing Bad Channels.');
%%% IF THERE ARE NO BAD CHANNELS TO REMOVE --------------------------------
    else
        fprintf('\n\nSUBJECT: %s\n has no bad channels. \n', subj_match.pid{ss});
        save([ssSubjPath(ss) '/ft_meg_fullyProcessed'],'data','-v7.3')
    end
    
%% turn off diary
right_now = clock;
fprintf('%d:%d:%02.f       Done running **%s**.\n', ...
    right_now(4:6), mfilename)
diary off
    
end