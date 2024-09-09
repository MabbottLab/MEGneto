function Step2_ICA(config, pid, run_check_or_fix, visitnum)

% Step2_ICA will:
%       1) Carry out ICA on preprocessed data ("run")
%       2) Open the interactive component checker GUI ("check"), or
%       3) Backproject components to be removed and run channel repair ("fix")
% 
% NOTES:
%   - Check participants who had excessive head motion or excessive numbers
%   of bad channels and remove them from this step if necessary. 
%   - Need to omit bad channels from ica - if not you will get complex 
%   numbers. Because during repair channels procedure bad channels are 
%   repaired according to neighbours, thus the new ones are not unique (no 
%   independent components).
%
% INPUTS:
%   > config: 
%       struct, configured in your "main" script with
%       all analysis parameters and options and paths
%   > pid:
%       string, participant ID used to build the paths
%       to relevant data files and output folders
%   > run_check_or_fix:
%       string, either "run", "check", or "fix" corresponding to which
%       utility you'd like to use
%   > visitnum:
%       int, visit number if longitudinal, optional arg)

% OUTPUTS:
%   > step2_data_fullyProcessed.mat: 
%       post-ICA input into step3
%   > step2_icaComponents.mat:
%       output of running ICA ("run") on step1 data output
%   > step2_badComp.csv:
%       CSV list of component #s to be rejected, outputted by the "check"
%       option
%
% See also: FT_DENOISE_SYNTHETIC, FT_RESAMPLEDATA, FT_SELECTDATA, FT_COMPONENTANALYSIS 

% Last updated by: Julie Tseng, 2020-01-08
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.
%% SETUP: LOAD THINGS

if exist('visitnum', 'var')
    this_output = [config.meta.project_path '/' config.meta.analysis_name '/' ...
                    pid '/' sprintf('ses-%.2d', visitnum)]; % indicate subject-specific output folder path
else
    this_output = [config.meta.project_path '/' config.meta.analysis_name '/' pid]; 
end

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
        writematrix(bad_comp, [this_output '/step2_badComp.csv']) 

%% regress ICA noise components and fix bad channels
    elseif run_check_or_fix == "fix"
        
        load([this_output '/step2_icaComponents.mat'])
        load([this_output '/step1_data_clean.mat'])
        load([this_output '/step2_badComp.csv'])
        
        % first, ICA component regress (if there are any)
        if ~isempty(step2_badComp)
            cfg = [];
            cfg.component = step2_badComp;
            data_clean = ft_rejectcomponent(cfg, comp, data_clean);
        end
        
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
