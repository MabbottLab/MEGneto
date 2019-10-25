function fcp_3_ChannelRepair(paths)

% After a human identifies ICA components that correspond to artifacts,
% this code will back-project leftover ICA components to data, as well as
% repair bad channels with the average neighbour signal. 

config = load_config(paths, paths.name);
config = config.config;
step = 'fcp3';
subj_ds = load_participants(paths,step);
pids = readtable(paths.all_subj_pids);
[subj_match, failure] = ds_pid_match(paths,step);
ssSubjPath = @(x) paths.(subj_match.pid{x});
rawdata = cellfun(@(x) [paths.rawdata '/' x], subj_match.ds, 'UniformOutput', false);

if isempty(subj_match)
    error('No participants selected')
end

write_match_if_not_empty(paths,step);

if length(unique(subj_match.pid)) ~= length(subj_match.pid)
    error('More than one ds per participant!')
end

fcp1_output = load_config(paths,'fcp1_output');
fcp2_output = loadjson([paths.anout_grp '\fcp2_5_output.json']);
fcp2_output = recursive_json_struct_string_to_func(fcp2_output);

rangeOFsubj = 1:height(subj_match);

disp('Starting channel repair...');
for ss = rangeOFsubj
    load([ssSubjPath(ss) '/' fcp2_output.preprocessedData_cfg],'-mat','data');

    % remove and repair bad channels
    if config.cleaningOptions.rmBadChannels == 1

        disp('Finding neighbours ...')

        cfg = [];
        cfg.method        = 'distance';%'distance', 'triangulation' or 'template'
        cfg.neighbourdist = 5;%number, maximum distance between neighbouring sensors (only for 'distance')
        cfg.template      = 'ctf151_neighb.mat'; % name of the template file, e.g. CTF275_neighb.mat
        %         cfg.layout        = filename of the layout, see FT_PREPARE_LAYOUT
        %         cfg.channel       = channels for which neighbours should be found
        %         cfg.feedback      = 'yes' or 'no' (default = 'no')
        neighbours = ft_prepare_neighbours(cfg, data);

        disp('Repairing bad channels ...')
        cfg = [];
        %         cfg.method         = 'weighted'; %'average', 'spline' or 'slap' (default = 'weighted')
        badChanns = fcp2_output.bad_chann(ss);
        badChanns = cellstr('M' + string(split(badChanns, 'M')));
        cfg.badchannel     = badChanns(2:length(badChanns));
        %cfg.missingchannel = cell-array, see FT_CHANNELSELECTION for details
        cfg.neighbours     = neighbours; %neighbourhood structure, see also FT_PREPARE_NEIGHBOURS
        %cfg.lambda         = regularisation parameter (default = 1e-5, not for method 'distance')
        %cfg.order          = order of the polynomial interpolation (default = 4, not for method 'distance')
        cfg.senstype = 'meg';
        data = ft_channelrepair(cfg, data);
        save([ssSubjPath(ss) '/ft_meg_fullyProcessed'],'data','-v7.3')

        disp('Done Removing Bad Channels.');

    end
end