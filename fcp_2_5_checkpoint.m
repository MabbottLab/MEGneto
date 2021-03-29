function fcp_2_5_checkpoint(paths)

% FCP_2_5_CHECKPOINT is an interactive piece that displays the timeseries
% associated with each ICA component, repeating for each participant. As
% data is displayed, the user will browse the components and then enter the
% ones that should be regressed from the data.
% 
% NOTES:
%   - This code will pull the subject csv from fcp_2 as all fcp_2
%   participants should have been analyzed. 
%
% INPUTS:
%   paths               =   struct defining paths to data, participant
%                           folders, analysis folders, config files, etc. 
%
% OUTPUTS:
%   fcp2_5_output       = struct with locations of output files
%     .data_noisecorr   = mat file with gradients accounted for
%     .preprocessedData_cfg = ft configuration
%     .ICAcomp_cfg      = JSON of ICA components
%     .data_icacomp     = mat file of ICA components
%     .ICA_badcomp      = record of user-specified bad components

%
% See also: DS_PID_MATCH, WRITE_MATCH_IF_NOT_EMPTY, FT_REJECTCOMPONENT,
% DISP_ICA_CHANS, FT_DATABROWSER

% Last updated by: Julie Tseng, 2020-01-08
%   This file is part of MEGneto, see https://github.com/SonyaBells/MEGneto
%   for the documentation and details.

%% SETUP

% logging file
right_now = clock;
log_filename = [paths.conf_dir '/log_' sprintf('%02.f%02.f%02.f', right_now(1:3))];
diary(log_filename)
fprintf('\n\n%02.f:%02.f:%02.f       Now running **%s**.\n', ...
    right_now(4:6), mfilename)

% load config JSON with analysis parameters
config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp2_5';

% check for matched MRI and MEG data
subj_match = ds_pid_match(paths,step);
ssSubjPath = @(x) paths.(subj_match.pid{x});
if isempty(subj_match) % if there are no full sets of data
    error('No participants selected')
end

% save the final participants list to the subj_match_fcp1 CSV
write_match_if_not_empty(paths,step);

% check for multiple *.ds folders for participants
if length(unique(subj_match.pid)) ~= length(subj_match.pid)
    error('More than one ds per participant!')
end

% load outputs from fcp_1 and fcp_2
fcp2_output = loadjson([paths.anout_grp '/fcp2_output.json']);
fcp2_output = recursive_json_struct_string_to_func(fcp2_output);

%% INTERACTIVE DISPLAY

% get the number of participants
rangeOFsubj = 1:height(subj_match);

for ss = rangeOFsubj
    
%%% LOAD AND DISPLAY DATA -------------------------------------------------
    load([ssSubjPath(ss) '/' fcp2_output.data_noisecorr],'-mat','data_noisecorr'); 
    disp_ica_chans(ss, ssSubjPath, config, fcp2_output);
    load([ssSubjPath(ss) '/' fcp2_output.data_icacomp])

    next    = false;
    gotin   = false;
    skip    = false;
    while ~next
        disp(subj_match.pid{ss}) % what participant?
        while ~gotin
%%% QUERY USER RESPONSE ---------------------------------------------------
            bad_comp = input('Enter the components to be removed (ex. [2,5,12]) ''skip'', or ''disp'': ');
        %%% EMPTY RESPONSE %%%
            if isempty(bad_comp)                                % if nothing was entered
                if strcmpi(input('Show again? (y/n)','s'),'Y')  % and user wants to see it again
                    gotin = false;                              % remain within the while loop
                else                                            % otherwise
                    gotin = true;                               % assume that there are no components to be removed
                    skip = false; 
                    next = 'Y';
                end
        %%% SKIP OR OTHER CHARACTER STRING ENTERED %%%
            elseif ~isnumeric(bad_comp)     % if non-numeric values were entered
                if ischar(bad_comp)         % and the values are characters
                    switch bad_comp
                        case 'skip'         % and the characters spell out 'skip'
                            next = true;
                            skip = true; 
                            break           % proceed to the next participant
                    end
                end
                disp_ica_chans(ss, ssSubjPath, config, fcp2_output); % DISPLAY
                gotin = false;
        %%% NUMERIC BAD COMPONENTS %%%
            else % when values are entered
                gotin   = true;
                % confirm bad components
                next    = upper(input(['Are these the bad components?: [' ...
                    num2str(bad_comp) ...
                    '] ([y]es/[n]o/[s]how components) '], 's'));
            end
        end
        gotin = false;
    %%% CONFIRM CHOSEN COMPONENTS %%%
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
%%% REMOVE BAD COMPONENTS -------------------------------------------------
    if ~skip
        % record bad components
        fcp2_output.bad_comp{ss,1} = {subj_match.pid{ss} bad_comp};

        % remove the bad components and backproject the data
        cfg             = [];
        cfg.component   = bad_comp; % to be removed component(s)
        data            = ft_rejectcomponent(cfg, comp, data_noisecorr);

        % save cleaned data
        save([ssSubjPath(ss) '/' fcp2_output.preprocessedData_cfg],'data','-v7.3');
        close all
    else % skip = 1, aka no bad components, re-save 
        fcp2_output.bad_comp{ss,1} = {subj_match.pid{ss} bad_comp}; % write "skip" in output to show no bad components were detected
        data            = data_noisecorr;
        save([ssSubjPath(ss) '/' fcp2_output.preprocessedData_cfg],'data','-v7.3');
        close all
    end

    % save a JSON copy of the components
    right_now = clock;
    save_to_json(bad_comp, sprintf('%s/ICA_badcomp_%02.f%02.f%02.f_%02.f%02.f.json', ...
        ssSubjPath(ss), right_now(1:5)));
end

%%% RECORD KEEPING --------------------------------------------------------

% save all bad ICA components
if config.cleaningOptions.artifact.icaClean == 1
    all_bad_comp = fcp2_output.bad_comp;
    save_to_json(all_bad_comp,[paths.anout_grp '/' fcp2_output.ICAcomp_cfg]);
end

% save fcp_2_5 output
disp('Saving fcp_2_5 output...');
save_to_json(fcp2_output, [paths.anout_grp '/fcp2_5_output.json'])
disp('Done.\n')

%% turn off logging

right_now = clock;
fprintf('%02.f:%02.f:%02.f       Done running **%s**.\n', ...
    right_now(4:6), mfilename)
diary off

end

%%% DISPLAY ICA CHANNEL MINI FUNCTION -------------------------------------
function disp_ica_chans(ss, ssSubjPath, config, fcp2_output)
    close all

    load([ssSubjPath(ss) '/' fcp2_output.data_icacomp],'-mat','comp');
    cfg             = [];
    cfg.channel     = [1:5]; % components to be plotted
    cfg.viewmode    = 'component';
    cfg.layout      = config.filteringParameters.CTFlayout;
    cfg.axisfontsize = 8;
    cfg.linewidth   = 0.2;
    cfg.plotlabels  = 'yes';
    cfg.position    = [300 200 1500 800];
    ft_databrowser(cfg, comp);
end