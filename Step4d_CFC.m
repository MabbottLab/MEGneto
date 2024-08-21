 function Step4d_CFC(config, pid, visitnum)
% CFC calculates the mean vector length (a phase amplitude coupling
% measure) between regions. 
%
% Configuration structure (config.step4d) has to contain:
%       lowFreqband: double with low frequency band (Hz) (ex: [4 8])
%       highFreqband: double with high frequency band (Hz) (ex: [8 12]) 
%
% Additional configuration options: 
%       chanlow: cell array with channel names to be
%       filtered into low frequency band (ex: {'chan1','chan2'} )
%       chanhigh:array with channel names to be
%       filtered into high frequency band
%
% Saves:
%       .mat file containing data (cfc values chanlow x chanhigh), labels
%       for chanlow, chanhigh
% 
% Notes:
%       - does not keep trials, averages them out 
%       - whenever a list of channels (chanlow or chanhigh) isn't specified, 
%        will default to using all the ROIs in the original data
%% SETUP
    if exist('visitnum', 'var')
        this_output = [config.meta.project_path '/' config.meta.analysis_name '/' pid '/' sprintf('ses-%.2d', visitnum)]; % indicate subject-specific output folder path
    else
        this_output = [config.meta.project_path '/' config.meta.analysis_name '/' pid]; 
    end

    load([this_output '/out_struct.mat'])
    load([this_output '/step3_data_roi.mat'])
%% CHECKS
    % checking for missing, required fields
    if ~isfield(config.step4d,'lowFreqband')|| ~isfield(config.step4d,'highFreqband')
        error('Please specify both high and low frequency band ranges');
    end
    
    % checking frequency channels 
    chan_fields = {'chanlow','chanhigh'};
    
    for i=1:length(chan_fields)
        if isfield(config.step4d,chan_fields(i)) 
            missing_chan = ~(ismember(config.step4d.(sprintf('%s',chan_fields{i})),data_roi.label)); % returns 1 for missing channels
            if sum(missing_chan)>0 %if there's at least 1 missing channel,display error with channel names
                missing_chan_names = strjoin(string(config.step4d.(sprintf('%s',chan_fields{i}))(missing_chan)));
                error(sprintf('The following channels for %s do not exist in data: %s.',chan_fields{i},missing_chan_names));
            end
            chans{i} = sort(ft_channelselection(config.step4d.(sprintf('%s',chan_fields{i})), data_roi.label));
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
    cfg          = [];
    cfg.bpfilter = 'yes'; 
    cfg.bpfreq   = config.step4d.lowFreqband;
    LFdata       = ft_preprocessing(cfg, data_roi);

    cfg          = [];
    cfg.bpfilter = 'yes'; 
    cfg.bpfreq   = config.step4d.highFreqband;
    HFdata       = ft_preprocessing(cfg, data_roi);
    
%% CONNECTIVITY CALCULATIONS
    ncomb = size(labelcmb,1); %number label combinations
    % saving names for plotting
    cfc.LFlabel = LFchan;
    cfc.HFlabel = HFchan; 

    % actual computation
    fprintf('CFC Connectivity Calculations');
    for i=1:nchanHF
        for j=1:nchanLF
            for k=1:ntrial
                % grab data for kth trial at current combination of
                % channels
                chandataLF = LFdata.trial{k}(strcmp(LFdata.label,LFchan(j)),:);
                chandataHF = HFdata.trial{k}(strcmp(HFdata.label,HFchan(i)),:);
                % abs = take magnitude of vector 
                % hilbert to get analytic signal
                cfc.data(j,i,k) = abs(mvl_calc(hilbert(chandataLF),hilbert(chandataHF)));
            end
        end
    end
    
    cfc.data = squeeze(mean(cfc.data,3));
    
    %%% saving data
    save([this_output sprintf('/step4d_CFC_%d_%d_vs_%d_%d.mat',config.step4d.lowF(1),config.step4d.lowF(2),config.step4d.highF(1),config.step4d.highF(2))], 'cfc'); 
end