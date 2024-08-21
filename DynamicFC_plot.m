function DynamicFC_plot(config)
% DynamicFC_plot creates ROI x ROI heatmap of dfc values for all specified
% participants and visits, from data calculated by Step4c_DynamicFC
%
% Configuration structure (config.DFCplot) must contain:
%       pids: patient ids to be included in heatmap
%       visitnum: cell array containing visit numbers for each participants
%       in order; participants do not have to have the same number of
%       visits
%           ex: config.DFCplot.visitnum = {[1] [1 2]} would use visit 1 for
%           pid 1 and visits 1 and 2 for pid 2
%       connmethod: connectivity method used in DFC calculations
%
% Saves:
%       dfc_avg.mat: average dfc values across all data (n_win x ROI x ROI
%       x freqband), pids and visitnums used to calculate values
%       .png: saves figures of heatmaps for all specified frequency bands
%       as avg_dfc_freq_%d_%d.png where %d is replaced by the bounds of the
%       frequency band 
% 
% Notes:
%       - works if given just one pid and one visit number  
%       - assumes all data contains same frequency bands, ROIs and windows  

    %%% CHECKS -----------------------------------------------------------
    % checking that required fields exist
    if ~isfield(config.DFCplot,'pids')
        error('Please specify patient ids');
    end
    
    if ~isfield(config.DFCplot,'visitnum')
        error('Please specify visit numbers');
    end
    
    if ~isfield(config.DFCplot,'connmethod')
        error('Please specify the connectivity method');
    end
    
    % checking that the number of pids correspond to the number of sets of
    % visitnumbers 
    if size(config.DFCplot.pids,2)~=size(config.DFCplot.visitnum,2)
        error('Please specify visit numbers for each participant');
    end
    
    %%% EXTRACTING DFC VALUES --------------------------------------------
    i=1; % dataset counter
    for p=1:length(config.DFCplot.pids)
        pid    = config.DFCplot.pids(p);
        visits = config.DFCplot.visitnum{p};
        
        for v=1:length(visits)
            % grab dfc matrix for that participant and visit number
            this_output = [config.meta.project_path '/' config.meta.analysis_name sprintf('/WMP_%.2d/ses-%.2d', pid, visits(v))]; 
            load([this_output '/step4c_connDFC.mat']);
            
            p_dfc = connDFC.(sprintf('%sspctrm',config.DFCplot.connmethod)); 
            
            % assuming all participants have the same number of
            % windows, ROIs and frequency bands analyzed, grab this info 
            if p==1 && v==1
                % get dimensions of matrix to store all dfc data
                dims    = size(p_dfc);
                n_win   = dims(1);
                n_ROI   = dims(2);
                n_freq  = dims(4);
                dfc_all = zeros(dims); % matrix to store all dfc values
                
                % extracting variables for plotting
                labels    = connDFC.label;
                freqbands = connDFC.freqbands;
                wins      = connDFC.win{1};
            end
            % n_win x ROI x ROI x freqband x data i 
            dfc_all(:,:,:,:,i) = p_dfc; 
            
            i = i+1;
        end
    end 
    
    % calculating average and saving
    dfc_avg.data     = squeeze(mean(dfc_all,5)); % average across all data
    dfc_avg.pids     = config.DFCplot.pids;
    dfc_avg.visitnum = config.DFCplot.visitnum;
    
    % save data to group folder if multiple pids; if just one participant
    % doesn't save (would be the same as step4c_connDFC.mat)
    if length(config.DFCplot.pids)>1
        % check if group folder to save group results exists, if not, create it
        group_folder = [config.meta.project_path '/' config.meta.analysis_name '/group_results/DynamicFC']; 
        
        if not(isfolder(group_folder))
            mkdir(group_folder);
        end
        
        save([group_folder '/dfc_avg.mat'], 'dfc_avg');
    end
            
    % -----------------------------------------------------------------------------
    % SAVING ROI X ROI MATRICES OVER TIME AS FIGURE FOR EACH
    % FREQUENCY BAND
    % -----------------------------------------------------------------------------     
    for fq = 1:size(freqbands,1) 
        % selecting average dfc for specific freq band
        dfc_avg_fq = dfc_avg(:,:,:,fq);
        
        % store max and min color values for uniform colorbar
        colorLims = [];

        % prevent figures from popping up 
        %set(0, 'DefaultFigureVisible','off'); 

        % move onto rows of 4 cols if more than 4 time windows
        if n_win > 4
            if mod(n_win,4)==0 % if n_win is divisible by 4
                n_row = n_win/4;
            else
                n_row = floor(n_win/4) + 1;
            end
            n_col = 4; 
        else
            n_row = 1;
            n_col = n_win;
        end 

        % create tiled layout figure for heatmaps
        fig = figure(fq); 
        tcl = tiledlayout(fig,n_row,n_col,'TileSpacing','Compact'); 
        h   = gobjects(n_row,n_col); 

        % plot heatmap for each window 
        for w = 1:n_win
            win_start = wins(w);
            win_end   = wins(w+1);

            % heatmap at window "w"
            nexttile(tcl);
            h(w) = heatmap(labels,labels,squeeze(dfc_avg_fq(w,:,:)));

            %%% formatting 
            if mod(w+3,4)~=0 % if the index+3 is not divisible by 4 (i.e. not in the first col)
                h(w).YDisplayLabels = nan(size(h(w).YDisplayLabels)); % remove ylabels
            end

            colorbar off
            h(w).NodeChildren(3).TickLabelInterpreter = 'none'; % prevent letter subscripts for ROI names
            set(h(w),'GridVisible','off','FontSize',6); 
            colorLims = [colorLims; h(w).ColorLimits];

            % label matrix with time window 
            title_text = compose([sprintf("%.3f to %.3f s",win_start,win_end)]); 
            h(w).Title = '\fontsize{10}' + title_text;
        end 
        
        % make colorbar constant for figure 
        globalColorLim = [min(colorLims(:,1)), max(colorLims(:,2))];
        set(h(1), 'ColorLimits', globalColorLim)
        
        ax             = axes(tcl,'visible','off','Colormap',h(1).Colormap,'CLim',globalColorLim);
        cb             = colorbar(ax);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
        cb.Layout.Tile = 'East';
        
        if length(config.DFCplot.pids)>1
            title(tcl,sprintf('Average Dynamic Functional Connectivity For Frequency %d-%d Hz',freqbands(fq,1),freqbands(fq,2)));
            saveas(fig,[group_folder sprintf('/avg_dfc_freq_%d_%d.png',freqbands(fq,1),freqbands(fq,2))])
        else
            title(tcl,sprintf('WMP_%.2d Dynamic Functional Connectivity For Frequency %d-%d Hz',pid,freqbands(fq,1),freqbands(fq,2)),'interpreter','none');
            saveas(fig,[this_output sprintf('/dfc_freq_%d_%d.png',freqbands(fq,1),freqbands(fq,2))])
        end
    end
end