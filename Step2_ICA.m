function Step2_ICA(config, pid, run_check_or_fix)

% same wrapper function for initial ICA component analysis, custom
% browsing, and regressing components

%% SETUP: LOAD THINGS

this_output = [config.meta.project_path '/' config.meta.analysis_name '/' pid];
load([this_output '/out_struct.mat']) % out variable

%% run ICA
    if run_check_or_fix == "run"
        if config.step2.icaClean == 1
    %%% ICA -------------------------------------------------------------------
            % first, handle bad channels if they exist
            % need to exclude the bad channels from ICA
            % load info about bad channels
            load([this_output '/step1_data_clean.mat']) % data_clean variable

            if ~isempty(out.step1.badChanDef.out)
                cfg         = [];
                cfg.channel = setdiff(data_clean.label, out.step1.badChanDef.out);         % specify included channels
                data_clean_rmBadCh = ft_selectdata(cfg, data_clean);     % select channel data

                % Run ICA
                cfg          = []; % set up config for ICA
                cfg.channel  = 'MEG';
                cfg.method   = 'fastica'; % default and uses the implementation from EEGLAB
                comp         = ft_componentanalysis(cfg, data_clean_rmBadCh); % run ICA
            else % if there are no bad channels, proceed with ICA
                cfg          = []; % set up config for ICA
                cfg.channel  = 'MEG';
                cfg.method   = 'fastica'; 
                comp         = ft_componentanalysis(cfg, data_clean); % run ICA
            end
            % save the ICA components
            save([this_output '/step2_icaComponents.mat'], 'comp')
        else
            disp("Tried to run ICA but you specified no in analysis config.")
        end
        
%% check through components interactively
    elseif run_check_or_fix == "check"
        load([this_output '/step2_icaComponents.mat'])
        
        close all
        cfg             = [];
        cfg.channel     = [1:5]; % components to be plotted
        cfg.viewmode    = 'component';
        cfg.layout      = 'CTF151.lay';
        cfg.axisfontsize = 8;
        cfg.linewidth   = 0.2;
        cfg.plotlabels  = 'yes';
        cfg.position    = [300 200 1500 800];
        ft_databrowser(cfg, comp);
        
        bad_comp = input('Enter the components to be removed in this format: [2, 5, 12]:');
        close all
        
        out.step2.ica_bad_comp = bad_comp;
        writematrix(bad_comp, [this_output '/step2_badComp.csv']) % note to JT: fix this

%% regress ICA noise components and fix bad channels
    elseif run_check_or_fix == "fix"
        
        load([this_output '/step2_icaComponents.mat'])
        load([this_output '/step1_data_clean.mat'])
        load([this_output '/step2_badComp.csv'])
        
        % first, ICA component regress (if there are any)
%        if ~isempty([this_output '/step2_badComp.mat'])
            cfg = [];
            cfg.component = step2_badComp;
            data_clean = ft_rejectcomponent(cfg, comp, data_clean);
%        end
        
        % then fix bad channels
        if ~isempty(out.step1.badChanDef.out)
            % prep neighbours
            cfg                 = [];
            cfg.method          = 'distance';
            cfg.neighbourdist   = 5;
            cfg.template        = 'ctf151_neighb.mat';
            neighbours = ft_prepare_neighbours(cfg, data_clean);
            
            % record some info
            out.step2.chanRepair.neighbourmethod = 'distance';
            out.step2.chanRepair.neighbourdist = 5;
            out.step2.chanRepair.template = 'ctf151_neighb.mat';
            out.step2.chanRepair.method = 'weighted';
            
            % actual channel repair w/ weighted avg of neighbours
            cfg                 = [];
            cfg.method          = 'weighted';
            cfg.badchannel      = out.step1.badChanDef.out;
            cfg.neighbours      = neighbours;
            cfg.senstype        = 'meg';
            data_clean          = ft_channelrepair(cfg, data_clean);
        end
        
        save([this_output '/step2_data_fullyProcessed'], 'data_clean')
        save([this_output '/out_struct.mat'], 'out')
    end

end
