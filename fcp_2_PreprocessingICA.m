function fcp_2_PreprocessingICA(paths)

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
fcp2_output.preprocessedData_grad =  'ft_meg_data_grad.json';
fcp2_output.ICAcomp_cfg           =  'ft_icacomp.json'; % save ica components
fcp2_output.subj_icafigure        =  @(ext) ['icafigure.', ext];
fcp2_output.subj_dir_rawchanfig   =  'raw_channel_figure';
fcp2_output.data_noisecorr        =  'data_noisecorr.mat';
fcp2_output.data_icacomp          =  'icacomponents.mat';
for pids = 1:length(subj_match.pid)
    mkdir([paths.(subj_match.pid{pids}) '/' fcp2_output.subj_dir_rawchanfig]);
end

% determine which subjects to process
rangeOFsubj = 1:height(subj_match);

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

fcp2_output.bad_chann = loadjson([paths.anout_grp '/group_rmBadChan.json']);
for ss = rangeOFsubj
    load([ssSubjPath(ss) '/' fcp2_output.data_noisecorr],'-mat','data_noisecorr');

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
        neighbours = ft_prepare_neighbours(cfg, data_noisecorr);

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
        data_noisecorr = ft_channelrepair(cfg, data_noisecorr);
        save([ssSubjPath(ss) '/' fcp2_output.data_noisecorr],'data_noisecorr','-v7.3')

        disp('Done Removing Bad Channels.');

    end

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

        % downsample the data to speed up the next step
        cfg = [];
        cfg.resamplefs = 300;
        cfg.detrend    = 'no';
        data_resampICA = ft_resampledata(cfg, data_noisecorr);

        % ICA - to find EOG and ECG
        % perform the independent component analysis (i.e., decompose the data)
        if  ~isempty(fcp2_output.bad_chann{ss})
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
            cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
            comp = ft_componentanalysis(cfg, data_resampICA);
        end
        save([ssSubjPath(ss) '/' fcp2_output.data_icacomp],'comp','-v7.3')


        % Identify artifacts
        % plot the components for visual inspection
        figure
        cfg = [];
        cfg.component = [1:20];       % specify the component(s) that should be plotted
        cfg.layout    = config.filteringParameters.CTFlayout; % specify the layout file that should be used for plotting
        cfg.comment   = 'no';
        ft_topoplotIC(cfg, comp)
        % edit display of components
        hf = gcf;
        hf.Position(1:2) = [10 10];
        hf.Position(3:4) = (800 / hf.Position(4)) .* hf.Position(3:4);
        %         disp('Close figure to move on');
        % save component as image
        print(hf, [ssSubjPath(ss) '/' fcp2_output.subj_icafigure('png')], '-dpng', '-r600');

        clear('comp')



    end

    close all
    fprintf('\nDone subject %s! \n',subj_match.pid{ss})

end

if interactive == true
    for ss = rangeOFsubj
        load([ssSubjPath(ss) '/' fcp2_output.data_noisecorr],'-mat','data_noisecorr');
        disp_ica_chans(ss, ssSubjPath, config, fcp2_output);
        load([ssSubjPath(ss) '/' fcp2_output.data_icacomp])
        %         pause(70) %pause for 4 minutes
        next = false;
        gotin = false;
        skip = false;
        while ~next
            disp(subj_match.pid{ss})
            while ~gotin
                bad_comp = input('Enter the components to be removed (ex. [2,5,12]) ''skip'', or ''disp'': ');
                if isempty(bad_comp)
                    if strcmpi(input('Show again? (y/n)','s'),'Y')
                        gotin = false;
                    else
                        gotin = true;
                        skip = false;
                        next = 'Y';
                    end
                elseif ~isnumeric(bad_comp)
                    if ischar(bad_comp)
                        switch bad_comp
                            case 'skip'
                                next = true;
                                skip = true; %Enter the string 'skip' to skip that participant
                                break
                        end
                    end
                    disp_ica_chans(ss, ssSubjPath, config, fcp2_output);
                    gotin = false;
                else
                    gotin = true;
                    next = upper(input(['Are these the bad components?: [' ...
                        num2str(bad_comp) ...
                        '] ([y]es/[n]o/[s]how components) '], 's'));
                end
            end
            gotin = false;
            if ~skip
                switch next
                    case 'Y'
                        next = true;
                    case 'N'
                        next = false;
                    case 'S'
                        next = false;
                        disp_ica_chans(ss, ssSubjPath, config, fcp2_output);
                    otherwise
                        next = false;
                        warning("Please enter y, n, or s after choosing components")
                end
            end
        end
        if ~skip
            % Get bad components from user
            fcp2_output.bad_comp{ss,1} = bad_comp;
            
            % Resample original dataset
            
            cfg = [];
            cfg.resamplefs = config.filteringParameters.sampleRate;
            cfg.detrend    = 'no';
            data_resamp = ft_resampledata(cfg, data_noisecorr);
            
            % remove the bad components and backproject the data
            cfg = [];
            cfg.component = bad_comp; % to be removed component(s)
            data = ft_rejectcomponent(cfg, comp, data_resamp);
            
            % save cleaned data
            save([ssSubjPath(ss) '/' fcp2_output.preprocessedData_cfg],'data','-v7.3');
            close all
        end
        
        save_to_json(bad_comp, [ssSubjPath(ss) '/ICA_badcomp.json'])
    end
    
    % save all bad components (just in case)
    if config.cleaningOptions.artifact.icaClean == 1
        all_bad_comp = fcp2_output.bad_comp;
        save_to_json(all_bad_comp,[paths.anout_grp '/' fcp2_output.ICAcomp_cfg]);
    end
    
    % save output
    disp('Saving...');
    save_to_json(fcp2_output,[paths.anout_grp '/fcp2_output.json'])
    disp('Done.')
    % % save p-structure
    % disp('Saving...');
    % save(p.paths.p_strct,'p','-mat','-v7');
    % disp('Done.');
    
end
end
function disp_ica_chans(ss, ssSubjPath, config, fcp2_output)
close all
imshow([ssSubjPath(ss) '/' fcp2_output.subj_icafigure('png')])
load([ssSubjPath(ss) '/' fcp2_output.data_icacomp],'-mat','comp');
cfg = [];
% cfg.blocksize = 1000;
cfg.channel = [1:5]; % components to be plotted
cfg.viewmode = 'component';
cfg.layout = config.filteringParameters.CTFlayout;
cfg.fontsize = 0.02;
cfg.axisfontsize = 8;
cfg.linewidth = 0.2;
ft_databrowser(cfg, comp);
% edit display
hf2 = gcf;
hf2.Position(1:2) = [300 200];
hf2.Position(3:4) = [1500,800];
end
