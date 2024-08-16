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
    
    % checking low frequency channels 
    if isfield(config.step4d,'chanlow')
        LFchan = sort(ft_channelselection(config.step4d.chanlow, data_roi.label));
        
        % checking if any of the inputted channels don't exist
        % if it doesn't, will get 0s 
        chanLF_not_exist = find(ismember(LFchan,data_roi.label)==0);
  
        % if there are channels that don't exist, looks for labels and
        % reports them 
        if ~isempty(chanLF_not_exist)
            chanLF_wrong = strjoin(string(config.step4d.chanlow(chanLF_not_exist)));
            error(sprintf('The following channels for low frequency do not exist in data: %s.',chanLF_wrong));
        end    
    % if chanlow not specified, automatically use all in data
    else 
        LFchan = sort(ft_channelselection('all', data_roi.label));
    end
    
    % repeat for high frequency channels 
    if isfield(config.step4d,'chanhigh')
        HFchan = sort(ft_channelselection(config.step4d.chanhigh, data_roi.label));
        chanHF_exist = find(ismember(HFchan,data_roi.label)==0); 
        if ~isempty(chanHF_exist)
            chanHF_wrong = strjoin(string(config.step4d.chanhigh(chanHF_exist)));
            error(sprintf('The following channels for high frequency do not exist in data: %s.',chanHF_wrong));
        end
    else
        HFchan = sort(ft_channelselection('all', data_roi.label));
    end
      
    ntrial = numel(data_roi.trial);
    nchanLF = numel(LFchan);
    nchanHF = numel(HFchan);

    % The line below puts lowF channels left col, highF on right, patterned
    % alternating lowF first; includes within channel 
    %{
    ex:
    chan1 chan1
    chan2 chan1
    chan1 chan2
    chan2 chan2
    %}
    labelcmb = ft_channelcombination({LFchan,HFchan},union(data_roi.label,data_roi.label),1,2);
    
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

    % create matrix that is lowchan x highchan x trials 
    cfc.data  = zeros(nchanLF,nchanHF,ntrial);     
    cfc.LFlabel = LFchan;
    cfc.HFlabel = HFchan; 
    % channel combination number 
    comb_num = 1;
    
    % actual computation
    fprintf('CFC Connectivity Calculations');
    for i=1:nchanHF
        for j=1:nchanLF
            for k=1:ntrial
                % grab data for kth trial at current combination of
                % channels
                chandataLF = LFdata.trial{k}(strcmp(LFdata.label,labelcmb{comb_num,1}),:);
                chandataHF = HFdata.trial{k}(strcmp(HFdata.label,labelcmb{comb_num,2}),:);
                % abs = take magnitude of vector 
                % hilbert to get analytic signal
                cfc.data(j,i,k) = abs(mvl_calc(hilbert(chandataLF),hilbert(chandataHF)));
            end
            comb_num = comb_num + 1;
        end
    end
    
    %%% plotting 
    cfc_plot = squeeze(mean(cfc.data,3));
    
    %%% saving data
    cfc.data = cfc_plot;
    save([this_output sprintf('/step4d_CFC_%d_%d_vs_%d_%d.mat',config.step4d.lowF(1),config.step4d.lowF(2),config.step4d.highF(1),config.step4d.highF(2))], 'cfc'); 
end