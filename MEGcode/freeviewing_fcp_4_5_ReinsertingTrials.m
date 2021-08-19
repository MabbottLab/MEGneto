%% Instructions
% run this first section of code (labelled section 1) with the CSV file 
% from step 4 excluding people with just one run of data (16, 17, P01). 
% The next section (labelled section 2) is for the remaining
% participants.

% This code will stack outputs of fcp4 on a per participant basis (i.e.,
% combine runs 1 and 2 for each participant that has 2 runs) and then
% insert blank trials for trials that were rejected in fcp1.

%%% ACTION NEEDED FROM THE USER
% This code uses the subj_fcp4.csv file to loop through participants.
% Please remove particiapnts 16, 17 and P01 from this csv file before
% running the code. You must also alter the paths in the first section of
% this code under the heading "Set some paths up".

%%% OUTPUTS
% Outputs from the fcp_4_5 function get saved in the run 1 analysis
% folders. From here on, please ensure the "paths" variable is set to
% be the one from the run 1 analysis, not run 2.
% The output is saved under participant specific folder under the name 
% 'revamped_beamforming_results.mat'.

%% Set some paths up (this should be altered depending on the user's path)
run1_analysis = '/home/dmatic/MEGProjects/analysis/OrderedClips_OIRMRun1_Jun23/analysis/'; % path to analysis folder for run 1 data
run2_analysis = '/home/dmatic/MEGProjects/analysis/OrderedClips_OIRMRun2_Jun25/analysis/'; % path to analysis folder for run 2 data

%%  Load config JSON with analysis parameters
config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp4'; % assumes that same ppts are used here as were used in the beamforming step

% check for matched MRI and MEG data
subj_match  = freeviewing_ds_pid_match(paths,step);

%% Remove ppts 16, 17, P01 from the subj_csv file


%% Check that paticipants have been selected
ssSubjPath  = @(x) paths.(subj_match.pid{x});
if isempty(subj_match) % if there are no full sets of data
    error('No participants selected')
end

%% SECTION 1
%% For each participant, load the missing trial indices and insert NaN rows
rangeOFsubj = 1:length(subj_match.ds);

for ss = rangeOFsubj % for each participant that has matched MEG/MRI data
    
    % load first half of fcp4 results (part 1 run 1)
    load([run1_analysis ssSubjPath(ss) '/atlas_beamforming_results_part1.mat']);
    disp('now loading part 1 run 1')
    disp(ssSubjPath(ss))
    
    catmatrix_part1 = catmatrix;
    coords_part1    = coords;
    srate_part1     = srate;
    
    % load second half of fcp4 results (part 2 run 1)
    % concatenate part 1 and part 2
    load([run1_analysis ssSubjPath(ss) '/atlas_beamforming_results_part2.mat']);
    disp('now loading part 2 run 1')
    disp(ssSubjPath(ss))
    
    catmatrix_run1 = [catmatrix_part1 catmatrix];
    coords_run1    = cat(1, coords_part1, coords);
    srate_run1     = srate;
    
    % fill in the missing trl indices for run 1
    trl_indices_1 = load([run1_analysis subj_match.pid{ss} '/rejected_trls.mat']);
    trl_indices_1 = trl_indices_1.trls_rejected;
    
    if ~isempty(trl_indices_1)
        for i = 1:length(trl_indices_1)
            current_col_1 = trl_indices_1(i); % where are we adding a new col?
            new_data_1 = NaN(size(catmatrix_run1(:,1,:)));
            catmatrix_run1 = [catmatrix_run1(:,1:current_col_1-1,:) new_data_1 catmatrix_run1(:,current_col_1:size(catmatrix_run1,2),:)];
        end
    end
    disp('dimension of catmatrix_run1')
    disp(size(catmatrix_run1,2))
    
    % load first half of fcp4 results (part 1 run 2)
    load([run2_analysis subj_match.pid{ss} '/atlas_beamforming_results_part1.mat']);
    disp('now loading part 1 run 2')
    disp(subj_match.pid{ss})
    
    catmatrix_part2 = catmatrix;
    coords_part2    = coords;
    srate_part2     = srate;
    
    % load second half of fcp4 results (part 2 run 2)
    % concatenate part 1 and part 2
    load([run2_analysis subj_match.pid{ss} '/atlas_beamforming_results_part2.mat']);
    disp('now loading part 2 run 2')
    disp(subj_match.pid{ss})
    
    catmatrix_run2 = [catmatrix_part2 catmatrix];
    coords_run2    = cat(1, coords_part2, coords);
    srate_run2     = srate;
    
    % fill in the missing trl indices for run 2
    trl_indices_2 = load([run2_analysis subj_match.pid{ss} '/rejected_trls.mat']);
    trl_indices_2 = trl_indices_2.trls_rejected;
    
    if ~isempty(trl_indices_2)
        for i = 1:length(trl_indices_2)
            current_col_2 = trl_indices_2(i); % where are we adding a new col?
            new_data_2 = NaN(size(catmatrix_run2(:,1,:)));
            catmatrix_run2 = [catmatrix_run2(:,1:current_col_2-1,:) new_data_2 catmatrix_run2(:,current_col_2:size(catmatrix_run2,2),:)];
        end
    end
    disp('dimensions of catmatrix_run2')
    disp(size(catmatrix_run2,2))
    
    % concatenate run 1 and run 2
    catmatrix = [catmatrix_run1 catmatrix_run2];
    
    % save srate
    if srate_run1 == srate_run2
        srate = srate_run1;
    else
        fprintf("mis-matched srates");
    end
    
    % stack and save coords
    coords = cat(1, coords_run1, coords_run2);% concatenate coords along dimension 1 
    
    % save it in the run 1 data folders
    save([run1_analysis subj_match.pid{ss} '/revamped_beamforming_results.mat'],'catmatrix', 'coords', 'srate', '-mat','-v7.3');
      
end

%% SECTION 2
%% complete the same process for the other participants who had only 1 run of data
outliars = ["OIRM16", "OIRM17", "OIRMP01"];

for i = 1:length(outliars)
    if outliars(i) == "OIRM16"
        data_part1 = load([run2_analysis '/OIRM16/atlas_beamforming_results_part1']);
        data_part2 = load([run2_analysis '/OIRM16/atlas_beamforming_results_part2']);
    
        % fill in the missing trl indices 
        trl_indices = load([run2_analysis '/OIRM16/rejected_trls.mat']);
        trl_indices = trl_indices.trls_rejected;
    
    elseif outliars(i) == "OIRM17"
        data_part1 = load([run1_analysis '/OIRM17/atlas_beamforming_results_part1']);
        data_part2 = load([run1_analysis '/OIRM17/atlas_beamforming_results_part2']);
    
        % fill in the missing trl indices 
        trl_indices = load([run1_analysis '/OIRM17/rejected_trls.mat']);
        trl_indices = trl_indices.trls_rejected;
    
    elseif outliars(i) == "OIRMP01"
        data_part1 = load([run2_analysis '/OIRMP01/atlas_beamforming_results_part1']);
        data_part2 = load([run2_analysis '/OIRMP01/atlas_beamforming_results_part2']);
    
        % fill in the missing trl indices 
        trl_indices = load([run2_analysis '/OIRMP01/rejected_trls.mat']);
        trl_indices = trl_indices.trls_rejected;
    end
    
    % data_part1
    catmatrix_part1 = data_part1.catmatrix;
    coords_part1    = data_part1.coords;
    srate_part1     = data_part1.srate;
    disp('catmatrix part 1')
    disp(size(catmatrix_part1,2))
    
    % data_part2
    catmatrix = [catmatrix_part1 data_part2.catmatrix];
    coords    = cat(1, coords_part1, data_part2.coords);
    srate_part2     = data_part2.srate;
    disp('catmatrix part 2')
    disp(size(data_part2.catmatrix,2))

    if ~isempty(trl_indices)
        for j = 1:length(trl_indices)
            current_col = trl_indices(j); % where are we adding a new col?
            new_data = NaN(size(catmatrix(:,1,:)));
            catmatrix = [catmatrix(:,1:current_col-1,:) new_data catmatrix(:,current_col:size(catmatrix, 2),:)];
        end
    end
    disp('catmatrix')
    disp(size(catmatrix,2))
    
    % add NaN ros to the end of ppt catmatrices that only have 1 run
    start_idx = size(catmatrix, 2) + 1;
    missing_length = 164 - size(catmatrix, 2);
    catmatrix(:,start_idx:164,:) = NaN;
    % save srate
    if srate_part1 == srate_part2
        srate = srate_part1;
    else
        fprintf("mis-matched srates");
    end
    
    % save it in the run 1 data folders
    if outliars(i) == "OIRM16"
        save([run1_analysis '/OIRM16/revamped_beamforming_results.mat'],'catmatrix', 'coords', 'srate', '-mat','-v7.3');
    elseif outliars(i) == "OIRM17"
        save([run1_analysis '/OIRM17/revamped_beamforming_results.mat'],'catmatrix', 'coords', 'srate', '-mat','-v7.3');
    elseif outliars(i) == "OIRMP01"
        save([run1_analysis '/OIRMP01/revamped_beamforming_results.mat'],'catmatrix', 'coords', 'srate', '-mat','-v7.3');
    end
    
end


