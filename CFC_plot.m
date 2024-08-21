function CFC_plot(config)
% CFC_plot creates ROI x ROI heatmap of CFC values for all specified
% participants and visits, from data calculated by Step4d_CFC
%
% Configuration structure (config.CFCplot) must contain:
%       pids: patient ids to be included in heatmap
%       visitnum: cell array containing visit numbers for each participants
%       in order; participants do not have to have the same number of
%       visits
%           ex: config.CFCplot.visitnum = {[1] [1 2]} would use visit 1 for
%           pid 1 and visits 1 and 2 for pid 2
%       lowFreqband: double with lower frequency band (in Hz)
%       highFreqband: double with higher frequency band (in Hz)
%
% Saves:
%       cfc_avg.mat: average cfc values across all data (ROI x ROI)
%       .png: saves figures of heatmaps 
% 
% Notes:
%       - works if given just one pid and one visit number  
%       - assumes all participants had same regions in low vs high
%       frequency 

    %%% CHECKS
    % checking that required fields exist
    if ~isfield(config.CFCplot,'pids')
        error('Please specify patient ids');
    end
    
    if ~isfield(config.CFCplot,'visitnum')
        error('Please specify visit numbers');
    end
    
    if ~isfield(config.CFCplot,'lowFreqband')
        error('Please specify the low frequency band');
    end
    
    if ~isfield(config.CFCplot,'highFreqband')
        error('Please specify the high frequency band');
    end
    
    % checking that the number of pids correspond to the number of sets of
    % visitnumbers 
    if size(config.CFCplot.pids,2)~=size(config.CFCplot.visitnum,2)
        error('Please specify visit numbers for each participant');
    end
    
    LF1= config.CFCplot.lowFreqband(1);
    LF2 = config.CFCplot.lowFreqband(2);
    HF1 = config.CFCplot.highFreqband(1);
    HF2 =config.CFCplot.highFreqband(2);

    % extracting CFC values for participants
    i=1; % dataset counter
    for p=1:length(config.CFCplot.pids)
        pid = config.CFCplot.pids(p);
        visits=config.CFCplot.visitnum{p};
        
        for v=1:length(visits)
            % grab cfc matrix for that participant and visit number
            this_output = [config.meta.project_path '/' config.meta.analysis_name sprintf('/WMP_%.2d/ses-%.2d', pid, visits(v))];
            load([this_output sprintf('/step4d_CFC_%d_%d_vs_%d_%d.mat',LF1,LF2,HF1,HF2)]);
            
            p_cfc = cfc.data;
            
            % this is assuming all participants have the ROIs in low and
            % high frequency
            if p==1 && v==1
                % get dimensions of matrix to store all cfc data
                nLF = numel(cfc.LFlabel);
                nHF = numel(cfc.HFlabel);
                cfc_all = zeros(nLF,nHF); % matrix to store all cfc values
                
                % extracting variables for plotting
                LFlabel = cfc.LFlabel;
                HFlabel = cfc.HFlabel;
            end
            % ROI x ROI x data i 
            cfc_all(:,:,i) = p_cfc;
            
            i = i+1;
        end
    end 
    
    % taking average across trials and saving
    cfc_avg = squeeze(mean(cfc_all,3));
    
    if length(config.CFCplot.pids)>1
        % check if group folder to save group results exists, if not, create it
        group_folder = [config.meta.project_path '/' config.meta.analysis_name '/group_results/CFC']; 
        
        if not(isfolder(group_folder))
            mkdir(group_folder);
        end
        
        save([group_folder '/cfc_avg.mat'], 'cfc_avg');
    end
    
    % -----------------------------------------------------------------------------
    % SAVING ROI X ROI MATRICES OVER TIME AS FIGURE 
    % -----------------------------------------------------------------------------
    fig = figure(1);
    h   = heatmap(LFlabel, HFlabel, cfc_avg');
    h.NodeChildren(3).TickLabelInterpreter = 'none'; % prevent letter subscripts for ROIs
    h.XLabel = 'Low Frequency Regions';
    h.YLabel = 'High Frequency Regions';
    
    if length(config.DFCplot.pids)>1
        h.Title = sprintf('Average CFC Between %d-%d Hz and %d-%d Hz',LF1,LF2,HF1,HF2);
        saveas(fig,[group_folder sprintf('/avg_CFC_%d_%d_vs_%d_%d.png',LF1,LF2,HF1,HF2)]);
    else
        h.Title = sprintf('WMP_%.2d CFC Between %d-%d Hz and %d-%d Hz',pid,LF1,LF2,HF1,HF2);
        saveas(fig,[this_output sprintf('/CFC_%d_%d_vs_%d_%d.png',LF1,LF2,HF1,HF2)]);
    end
end
