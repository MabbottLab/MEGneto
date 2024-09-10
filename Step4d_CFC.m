 function Step4d_CFC(config, pid, visitnum)
%
% Step4d_CFC is a cross-frequency coupling function that calculates the 
% mean vector length (a phase amplitude coupling measure) between regions. 
%
% INPUTS:
%   > config:
%       struct, configured in your "main" script with all analysis parameters 
%       and paths. Specifically, need the following fields:
%           > config.step4d.lowFreqband: double with low frequency band (Hz) (ex: [4 8])
%           > config.step4d.highFreqband: double with high frequency band (Hz) (ex: [8 12]) 
%
% Additional configuration options: 
%       > config.step4d.chanlow: 
%           cell array with channel names to be filtered into low frequency band (ex: {'chan1','chan2'} )
%       > config.step4d.chanhigh:
%           array with channel names to be filtered into high frequency band
%
% OUTPUTS:
%       > step4d_CFC_lowfreq_vs_highfreq.mat file containing data (cfc values chanlow x chanhigh), labels
%       for chanlow, chanhigh
% 
% Notes:
%       - does not keep trials, averages them out 
%       - whenever a list of channels (chanlow or chanhigh) isn't specified, 
%        will default to using all the ROIs in the original data
%
% Last updated by: Bianca Ha and Julie Tseng, 2024-09-09
%   This file is part of MEGneto, see https://github.com/MabbottLab/MEGneto
%   for the documentation and details.
%
%% SETUP
    if exist('visitnum', 'var')
        this_output = [config.meta.project_path '/' config.meta.analysis_name '/' ...
            pid '/' sprintf('ses-%.2d', visitnum)]; % indicate subject-specific output folder path
    else
        this_output = [config.meta.project_path '/' config.meta.analysis_name '/' pid]; 
    end

    load([this_output '/out_struct.mat'])
    load([this_output '/step3_data_roi.mat'])
    
    if isstring(data_roi.label)
        data_roi.label = cellstr(data_roi.label);
    end
    
%% CHECK FOR PARAMETER EXISTENCE IN CONFIG

    % did they specify the low/high frequency band ranges?
    if ~isfield(config.step4d,'lowFreqband')|| ~isfield(config.step4d,'highFreqband')
        error('Please specify both high and low frequency band ranges');
    end
    
    % checking frequency channels 
    chan_fields = {'chanlow','chanhigh'};
    
    for i = 1:length(chan_fields) % for chanlow and chanhigh
        if isfield(config.step4d, chan_fields(i)) % check for whether the config.step4d has selected sources
            % check whether the selected source is available in the
            % data_roi structure; returns 1 if missing
            missing_chan = ~(ismember(config.step4d.(sprintf('%s',chan_fields{i})),data_roi.label));
            
            % handle missing specified sources/channels
            if sum(missing_chan) > 0
                missing_chan_names = strjoin(string(config.step4d.(sprintf('%s',chan_fields{i}))(missing_chan)));
                error(sprintf('The following channels for %s do not exist in data: %s.',chan_fields{i},missing_chan_names));
            else % else, everything is fine and can proceed
                chans{i} = sort(ft_channelselection(config.step4d.(sprintf('%s',chan_fields{i})), data_roi.label));
            end
        else %if no channels specified, use all ROIs by default
            chans{i} = sort(ft_channelselection('all', data_roi.label));
        end
    end
    
    %grab ROIs 
    LFchan = chans{1};
    HFchan = chans{2};
      
    ntrial  = numel(data_roi.trial);
    nchanLF = numel(LFchan);
    nchanHF = numel(HFchan);
    
%% PREPROCESS DATA
    % low frequency band data filtering
    cfg          = [];
    cfg.bpfilter = 'yes'; 
    cfg.bpfreq   = config.step4d.lowFreqband;
    LFdata       = ft_preprocessing(cfg, data_roi);

    % high frequency band data filtering
    cfg          = [];
    cfg.bpfilter = 'yes'; 
    cfg.bpfreq   = config.step4d.highFreqband;
    HFdata       = ft_preprocessing(cfg, data_roi);
    
%% CONNECTIVITY CALCULATIONS
    % saving names for plotting
    cfc.LFlabel = LFchan;
    cfc.HFlabel = HFchan; 

    % actual computation
    fprintf('CFC Connectivity Calculations');
    for i = 1:nchanHF % for each high frequency channel/source
        for j = 1:nchanLF % and for each low frequency channel/source
            for k = 1:ntrial % and for each trial
                % grab data for kth trial at current combination of channels
                chandataLF = LFdata.trial{k}(strcmp(LFdata.label,LFchan(j)),:);
                chandataHF = HFdata.trial{k}(strcmp(HFdata.label,HFchan(i)),:);
                % abs = take magnitude of vector 
                % then hilbert(...) to get analytic signal
                cfc.data(j,i,k) = abs(mvl_calc(hilbert(chandataLF), ...
                                               hilbert(chandataHF)));
            end
        end
    end
    
    % take the mean across 
    cfc.data = squeeze(mean(cfc.data,3));
    
    % saving data
    save([this_output sprintf('/step4d_CFC_%d_%d_vs_%d_%d.mat',config.step4d.lowF(1),config.step4d.lowF(2),config.step4d.highF(1),config.step4d.highF(2))], 'cfc'); 
end