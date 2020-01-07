function fcp_2_PreprocessingICA(paths, included, skip_preprocessing)

%%% INPUTS %%% (update anonymous function's workspace if any edits to subject paths made)
%p.paths.epoching_info =     @(id,ext) [p.paths.output_dir(id),'epoching_',num2str(p.thr),'mmHM.',ext];

config = load_config(paths, paths.name);
config = config.config;
step = 'fcp2';
subj_ds = load_participants(paths,step);
pids = readtable(paths.all_subj_pids);
[subj_match, failure] = ds_pid_match(paths,step);
ssSubjPath = @(x) paths.(subj_match.pid{x});
% This will be modified when we settle on folder structure
% For now, temporary fix
% rawdata = cellfun(@(x) [paths.rawdata '/' x], subj_match.ds, 'UniformOutput', false);
rawdata = subj_match.ds;
if isempty(subj_match)
    error('No participants selected')
end

write_match_if_not_empty(paths,step);

if length(unique(subj_match.pid)) ~= length(subj_match.pid)
    error('More than one ds per participant!')
end

fcp1_output = load_config(paths,'fcp1_output');

%%% OUTPUT %%%
fcp2_output.preprocessedData_cfg  =  'ft_meg_data_cfg.mat';
fcp2_output.ICAcomp_cfg           =  'ft_icacomp.json'; % save ica components
fcp2_output.data_noisecorr        =  'data_noisecorr.mat';
fcp2_output.data_icacomp          =  'icacomponents.mat';

% determine which subjects to process
rangeOFsubj = 1:height(subj_match);

if ~skip_preprocessing
    disp('Starting Preprocessing...');
    for ss = rangeOFsubj %1:rangeOFsubj

        fprintf('\n\n==================================\nSUBJECT: %s\n', subj_match.pid{ss});

        % check for epoching file
        if exist([ssSubjPath(ss) '/' fcp1_output.trial_cfg], 'file')

            if config.task.isRest == false
                samples = loadjson([ssSubjPath(ss) '/' fcp1_output.trial_cfg]);
                cfg=[];
                cfg.trl = samples.trl;
                disp('Epoched file found');
                disp([ssSubjPath(ss) '/' fcp1_output.trial_cfg]);
            else
                if config.cleaningOptions.artifact.rmNoisyTrls == 0
                    samples = loadjson([ssSubjPath '/' fcp1_output.trial_cfgHM]);
                    cfg=[];
                    cfg.trl = samples.trl;
                    disp('Epoched file found');
                    disp([ssSubjPath(ss) '/' fcp1_output.trial_cfgHM]);
                elseif config.cleaningOptions.artifact.rmNoisyTrls == 1
                    samples = loadjson([ssSubjPath(ss) '/' fcp1_output.trial_cfg]);
                    cfg=[];
                    cfg.trl = samples.trl;
                    disp('Epoched file found');
                    disp([ssSubjPath(ss) '/' fcp1_output.trial_cfg]);
                end
            end
        else
            error('Epoched mat %s was not found', [ssSubjPath(ss) '/' fcp1_output.trial_cfg]);
        end

        % define data to import and filter
        cfg.dataset = [paths.rawdata '/' subj_match.ds{ss}];

        disp('Filtering Data....');

        cfg.channel = config.filteringParameters.channel;
        cfg.dftfilter = config.filteringParameters.dftfilter;
        cfg.dftfreq = config.filteringParameters.dftfreq;
        cfg.bpfilter = config.filteringParameters.bpfilter;
        cfg.bpfreq = config.filteringParameters.bpfreq;
        cfg.bpfiltord = config.filteringParameters.bpfiltord;
        cfg.continuous         = 'yes';
        % load data!
        dataFiltered = ft_preprocessing(cfg);

        % if modified gradiometer positions exist, put them in
        if exist([ssSubjPath(ss) '/' fcp1_output.grad_cfg] , 'file')
            fprintf('Loading gradiometer configuration from file...\n');
            dataFiltered.grad = loadjson([ssSubjPath(ss) '/' fcp1_output.grad_cfg]);
            dataFiltered.hdr.grad = dataFiltered.grad;
        else
            warning('No modified gradiometer definitions found!');
        end

        % synthetic 3rd order grads - for noise reduction in CTF no ref channels stored
        cfg = [];
        cfg.gradient = 'G3BR';
        data_noisecorr = ft_denoise_synthetic(cfg,dataFiltered);
        save([ssSubjPath(ss) '/' fcp2_output.data_noisecorr],'data_noisecorr','-v7.3')

    end
end

fcp2_output.bad_chann = loadjson([paths.anout_grp '/group_rmBadChan.json']);


for ss = rangeOFsubj

    if config.cleaningOptions.artifact.icaClean == 0

        % Resample all datasets
        cfg = [];
        cfg.resamplefs = config.filteringParameters.sampleRate;
        cfg.detrend    = 'no';
        data_resamp = ft_resampledata(cfg, data_noisecorr);

        % save Data file
        data = data_resamp;
        save([ssSubjPath(ss) '/' fcp2_output.preprocessedData_cfg],'data','-v7.3');

    elseif config.cleaningOptions.artifact.icaClean == 1
        
        disp('Running ICA');
        
        if skip_preprocessing
            load([ssSubjPath(ss) '/' fcp2_output.data_noisecorr])
        end
        
        % downsample the data to speed up the next step
        cfg = [];
        cfg.resamplefs = 300;
        cfg.detrend    = 'no';
        data_resampICA = ft_resampledata(cfg, data_noisecorr);

        % ICA - to find EOG and ECG
        % perform the independent component analysis (i.e., decompose the data)
        if  ~isempty(fcp2_output.bad_chann{included(ss)}{1})
            % get list of bad channels from fcp_1
            badChann = fcp2_output.bad_chann{included(ss)};
            badChanns = cellstr('M' + split(string(badChann), 'M'));
            %need to remove bad channels from ica - if not you will get
            %complex numbers. Because during repair channels procedure bad
            %channels are repaired according to neighbours, thus the new
            %ones are not unique (no independent components)
            % rmCh = strcat('-',fcp2_output.bad_chann{ss});
            % selchan = ft_channelselection({'all' rmCh{1,:}}, data_noisecorr.label);
            keepsies = setdiff(data_noisecorr.label, badChanns);
            cfg = [];
            cfg.channel = keepsies;
            data_resampICArmBadCh = ft_selectdata(cfg, data_resampICA);

            cfg = [];
            cfg.channel  = 'MEG';
            cfg.method = 'fastica'; % this is the default and uses the implementation from EEGLAB
            comp = ft_componentanalysis(cfg, data_resampICArmBadCh);

        else
            cfg = [];
            cfg.channel  = 'MEG';
            cfg.method = 'fastica'; % this is the default and uses the implementation from EEGLAB
            comp = ft_componentanalysis(cfg, data_resampICA);
        end
        save([ssSubjPath(ss) '/' fcp2_output.data_icacomp],'comp','-v7.3')

    end

    close all
    fprintf('\nDone subject %s! \n',subj_match.pid{ss})

end

disp('Saving...');
save_to_json(fcp2_output,[paths.anout_grp '/fcp2_output.json'])
disp('Done.')

end
