function fcp_2_5_checkpoint(paths)

% Interactive display selection of components that are noise

config = load_config(paths, paths.name);
config = config.config;
step = 'fcp2';
subj_ds = load_participants(paths,step);
pids = readtable(paths.all_subj_pids);
[subj_match, failure] = ds_pid_match(paths,step);
ssSubjPath = @(x) paths.(subj_match.pid{x});
% This will be modified when we settle on folder structure
% For now, temporary fix
rawdata = cellfun(@(x) [paths.rawdata '/' x], subj_match.ds, 'UniformOutput', false);

if isempty(subj_match)
    error('No participants selected')
end

write_match_if_not_empty(paths,step);

if length(unique(subj_match.pid)) ~= length(subj_match.pid)
    error('More than one ds per participant!')
end

fcp1_output = load_config(paths,'fcp1_output');
fcp2_output = loadjson([paths.anout_grp '\fcp2_output.json']);
fcp2_output = recursive_json_struct_string_to_func(fcp2_output);

rangeOFsubj = 1:height(subj_match);

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

%         % Resample original dataset
% 
%         cfg = [];
%         cfg.resamplefs = config.filteringParameters.sampleRate;
%         cfg.detrend    = 'no';
%         data_resamp = ft_resampledata(cfg, data_noisecorr);

        % remove the bad components and backproject the data
        cfg = [];
        cfg.component = bad_comp; % to be removed component(s)
        data = ft_rejectcomponent(cfg, comp, data_noisecorr);

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
save_to_json(fcp2_output,[paths.anout_grp '/fcp2_5_output.json'])
disp('Done.')

end

function disp_ica_chans(ss, ssSubjPath, config, fcp2_output)
    close all

    load([ssSubjPath(ss) '/' fcp2_output.data_icacomp],'-mat','comp');
    cfg = [];
    % cfg.blocksize = 1000;
    cfg.channel = [1:5]; % components to be plotted
    cfg.viewmode = 'component';
    cfg.layout = config.filteringParameters.CTFlayout;
    cfg.fontsize = 0.02;
    cfg.axisfontsize = 8;
    cfg.linewidth = 0.2;
    cfg.plotlabels = 'yes';
    ft_databrowser(cfg, comp);
    % edit display
    hf2 = gcf;
    hf2.Position(1:2) = [300 200];
    hf2.Position(3:4) = [1500,800];
end