function inspecting_results(paths, name, type)

% type.viz == 1 for imagesc plots by AAL region
% type.viz == 2 for plot of relative power over time in specific band
    % need also:
    % type.ppts = 'single', 'by_group', 'all'
        % if by_group chosen, need 
        % type.des_mat = matrix of participant breakdown
    % type.freq = [8 12], [13 29], [30 59], [60 100] (freq band)
    % type.freq_name = 'alpha, 'beta', ...
    % type.seeds = [1 2 19 20] (aal regions)
    % type.group = {[1, 2, 3, 4], [5, 6, 7]} (cell array of participant
    % groups with one cell per group
% type.viz == 3 for getting max power
    % type.show = true or false
    
config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp5';

% check for matched MRI and MEG data
subj_match  = ds_pid_match(paths,step);

% build list of indices assoc. w/ each group
groups = readtable([paths.conf_dir '/ParticipantCategories.xlsx']);
patients = ismember(subj_match.pid, groups.PATIENT);
ctrls = ismember(subj_match.pid, groups.TDC);

% timewindows
timewindows = -1.5:0.05:1.5;  

bands = cellstr(["alpha", "beta", "lgamma", "hgamma"]);
bands(2,:) = {[8 12], [13 29], [30 59], [60 100]};
    
load([paths.anout_grp '/' name]); % participant x AAL region x frequencies x timewindows
atlas = importdata('/mnt/sda/juanita/fieldtrip/template/atlas/aal/ROI_MNI_V4.txt');

if type.viz == 1 % show imagesc plots for each AAL region
    disp('Press any key to show next AAL region.')
    for region = 1:90
        imagesc(squeeze(mean(pow_spctrm(:,region,:,:),1)))
        ax = gca;
        ax.YTickLabel = (5:5:50).*2;
        ax.XTick = 11:10:61;
        ax.XTickLabel = timewindows(11:10:61);
        xlabel('Timewindows')
        ylabel('Frequencies 1-100Hz')
        title(sprintf('%s power spectrum', atlas.textdata{region,2}),'Interpreter','none');
        pause;
    end
elseif type.viz == 2 % plot
    
    % find frequencies within each band and that we also have calculated
    % power spectrum for
    if ~mod(type.freq(1),2) 
        these_freqs = (type.freq(1):2:type.freq(2))/2;
    else
        these_freqs = ((type.freq(1)+1):2:type.freq(2))/2;
    end
    
    % collapse across band
    % this_spctrm = max(pow_spctrm(:,:,these_freqs,:),[],3); % MAX
    this_spctrm = squeeze(nanmean(pow_spctrm(:,:,these_freqs,:), 3));    % MEAN
    
    if strcmp(type.ppts, 'all') % if displaying average across all
        this_spctrm = nanmean(this_spctrm, 1);
    elseif strcmp(type.ppts, 'by_group')
        this_spctrm = cat(1, ...
                            nanmean(this_spctrm(patients,:,:),1), ...
                            nanmean(this_spctrm(ctrls,:,:),1));
    end
        
    for seed = type.seeds
        for p = 1:size(this_spctrm,1)
            plot(squeeze(this_spctrm(p,seed,:)))
            ax = gca;
            ax.XTick = 11:10:61;
            ax.XTickLabel = timewindows(11:10:61);
            xlabel('Timewindows')
            ylabel('Relative power')
            if size(this_spctrm,1) == 1
                title(sprintf('%s power over time in %s, grouping: %s', ...
                    type.freq_name, atlas.textdata{seed,2}, type.ppts), ...
                    'Interpreter', 'none')
            elseif size(this_spctrm,1) == 2
                title(sprintf('%s power over time in %s, grouping: %s', ...
                    type.freq_name, atlas.textdata{seed,2}, groups.Properties.VariableNames{p}), ...
                    'Interpreter', 'none')
            else
                title(sprintf('%s power over time in %s, Participant %s', ...
                    type.freq_name, atlas.textdata{seed,2}, subj_match.pid{p}), ...
                    'Interpreter', 'none')
            end
            pause;
        end
    end
elseif type.viz == 3 % getting peak power etc.
    
    timewindows = -1.5:0.05:1.5;  
    all_max_tw = NaN(height(subj_match), size(pow_spctrm,2), length(bands));
    for f = 1:length(bands)
        % find frequencies within each band and that we also have calculated
        % power spectrum for
        if ~mod(bands{2,f}(1),2) 
            these_freqs = (bands{2,f}(1):2:bands{2,f}(2))/2;
        else
            these_freqs = ((bands{2,f}(1)+1):2:bands{2,f}(2))/2;
        end

        % get max value 
        this_spctrm = nanmax(pow_spctrm(:,:,these_freqs,:),[],3); % across freq band
        [~, max_tw] = nanmax(this_spctrm, [], 4); % across time
        
        all_max_tw(:,:,f) = max_tw;
        
        if type.show == true
                for region = 1:90
                    peak_power = zeros(height(subj_match), length(timewindows));
                    peak_power(sub2ind([52,61], [1:52]', max_tw(:,region))) = 1;
                    peak_power(:,31) = 0.5;
                    imagesc(peak_power)
                    ax = gca;
                    ax.XTick = 11:10:61;
                    ax.XTickLabel = timewindows(11:10:61);
                    xlabel('Timewindows')
                    ylabel('Participants')
                    title(sprintf('%s peak timewindow in %s band', atlas.textdata{region,2}, bands{1,f}),'Interpreter','none');
                    pause;
                end
        end
    end
    save([paths.anout_grp '/freqband_peak_timewindow.mat'],'all_max_tw','-mat','-v7.3')
    % output: participants x region x frequency x time
elseif type.viz == 4

    load([paths.anout_grp '/freqband_peak_timewindow.mat']); 

    for region = 1:90
        peak_power = zeros(height(subj_match), length(timewindows));
        peak_power(sub2ind([52,61], [1:52]', squeeze(all_max_tw(:,region, type.freq_band)))) = 1;
        peak_power(:,31) = 0.5;
        imagesc(peak_power)
        ax = gca;
        ax.XTick = 11:10:61;
        ax.XTickLabel = timewindows(11:10:61);
        xlabel('Timewindows')
        ylabel('Participants')
        title(sprintf('%s peak timewindow in %s band', atlas.textdata{region,2}, bands{1,type.freq_band}),'Interpreter','none');
        pause;
    end
elseif type.viz == 5
    diff_powspctrm = squeeze(...
                     nanmean(pow_spctrm(ctrls,:,:,:),1) - ...
                     nanmean(pow_spctrm(patients,:,:,:),1));
                 
    for region = 1:90
        imagesc(squeeze(diff_powspctrm(region,:,:)))
        ax = gca;
        ax.YTickLabel = (5:5:50).*2;
        ax.XTick = 11:10:61;
        ax.XTickLabel = timewindows(11:10:61);
        xlabel('Timewindows')
        ylabel('Frequencies 1-100Hz')
        title(sprintf('%s power spectrum, Controls-Patients', atlas.textdata{region,2}),'Interpreter','none');
        pause;
    end
end
    
    
    