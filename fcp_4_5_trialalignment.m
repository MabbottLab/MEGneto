function fcp_4_5_trialalignment(paths, trials_one, trials_two)

% FCP_4_5_TRIALALIGNMENT aligns participant trials that were conducted in
% two separate recordings. This also ensures that trials that were excluded
% from analysis (due to bad channels, excessive noise, etc) do not change
% the final order of participant trials.

%% SETUP
% load rejected trial numbers
load('/home/ckim/MEG_processing/analysis/trl_rejected.mat');
% load participants used for fcp_4
step        = 'fcp4';
subj_match  = ds_pid_match(paths,step);
ssSubjPath  = @(x) paths.(subj_match.pid{x});
rangeOFsubj = 1:length(subj_match.ds);
%% 	LOAD CAT MATRICES --------------------------------------------------
for ss = rangeOFsubj
    right_now = clock;
    fprintf('%02.f:%02.f:%02.f      Working on subject %s!\n', ...
        right_now(4:6), trl_rejected.Properties.RowNames{ss})
    
    % load first half of CATmatrix from run 1
    load([ssSubjPath(ss) '/atlas_beamforming_results1.mat']);
    
    catmatrix_run1  = catmatrix;
    coords_run1     = coords;
    srate_run1      = srate;
    
    % load second half of CATmatrix from run 1 and then stack the two 
    load([ssSubjPath(ss) '/atlas_beamforming_results2.mat']);
    
    catmatrix_run1  = [catmatrix_run1 catmatrix];
    coords_run1     = [coords_run1 coords];
    srate_run1      = srate;
    %% OPTIONAL: only if participants had two runs that need to be collated together
    if trials_two ~= 0
        %     load CATmatrix from run 2, excluding ppts 8, 16, and 17 since they only
        %     have one valid run
        % REMEMBER TO CHANGE PATH IF LOADING DATA FROM OTHER DATA FOLDERS
        if ss ~= 6 && ss~= 11 && ss~= 12
            load(['/home/ckim/MEG_processing/analysis/OIRMProcessingRun2/analysis/' trl_rejected.Properties.RowNames{ss} ...
                '/atlas_beamforming_results1.mat']);
            catmatrix_run2  = catmatrix;
            coords_run2     = coords;
            srate_run2      = srate;
            
            load(['/home/ckim/MEG_processing/analysis/OIRMProcessingRun2/analysis/' trl_rejected.Properties.RowNames{ss} ...
                '/atlas_beamforming_results2.mat']);
            catmatrix_run2  = [catmatrix_run2 catmatrix];
            coords_run2     = [coords_run2 coords];
            srate_run2      = srate;
        end
    end
%% Stacking cat matrices without removed trials
    catmatrix = catmatrix_run1;
    %% ONLY run if there were two runs and the optional section above was
    if trials_two ~= 0
        if ss ~= 6 && ss ~= 11 && ss~= 12
            if length(catmatrix) < length(catmatrix_run2)
                catmatrix(length(catmatrix)+1:length(catmatrix_run2), :, :) = NaN;
            elseif length(catmatrix) > length(catmatrix_run2)
                catmatrix_run2(length(catmatrix_run2)+1:length(catmatrix), :, :) = NaN;
            end
            catmatrix = [catmatrix catmatrix_run2];
        end
    end
    newcatmatrix = catmatrix;
    save([ssSubjPath(ss) '/compiledCATmatrix.mat'],'newcatmatrix', 'srate', '-mat','-v7.3')
end

%% Social content vector aligned with missing trials
% Social content vector associating 1s with social trials and 0s with
% non-social trials was made. Resultingly, a 164x1 column vector aligned
% with the content of each trial used.
% This portion ensures that for each participant, a new social content
% vector is made to align with removed trials due to signal artifacts, etc.

load([paths.anout '/conditionContentVector.mat']); % array named socialContent
content = conditionContent;                  
for ss = rangeOFsubj
    % ppt 6 had more trials in its first run than other ppts, which is why it is
    % under an elseif statement
    right_now = clock;
    fprintf('%02.f:%02.f:%02.f      Working on subject %s!\n', ...
        right_now(4:6), trl_rejected.Properties.RowNames{ss})
    if ~isempty(trl_rejected.Run1(ss)) && ss ~=6
        to_reject = cell2mat(trl_rejected.Run1(ss));
        tmp1 = content(1:trials_one, 1);
        for dd = 1:length(to_reject)
            tmp1(to_reject(dd)) = NaN;
        end
        remove_ind = isnan(tmp1);
        tmp1(remove_ind) = [];
    elseif ss == 6
        to_reject = cell2mat(trl_rejected.Run1(ss));
        tmp1 = content(1:131, 1);
        for dd = 1:length(to_reject)
            tmp1(to_reject(dd)) = NaN;
        end
        remove_ind = isnan(tmp1);
        tmp1(remove_ind) = [];
    else
        tmp1 = content(1:trials_one, 1);
    end
    if trials_two ~= 0
        if ~isempty(trl_rejected.Run2(ss)) && ss ~= 6 && ss~= 11 && ss~= 12
            to_reject = cell2mat(trl_rejected.Run2(ss));
            tmp2 = content(trials_one + 1:trials_one + trials_two, 1);
            for dd = 1:length(to_reject)
                tmp2(to_reject(dd)) = NaN;
            end
            remove_ind = isnan(tmp2);
            tmp2(remove_ind) = [];
            interest = [tmp1; tmp2];
        elseif ss == 6 || ss == 11 || ss == 12
            interest = tmp1;
        else
            tmp2  = content(trials_one + 1: trials_one+trials_two, 1);
            interest = [tmp1; tmp2];
        end
    end
    save([ssSubjPath(ss) '/conditionVector.mat'], 'interest')
end

right_now = clock;
fprintf('%02.f:%02.f:%02.f ============== Finished Processing ====================\n', ...
    right_now(4:6))
end
