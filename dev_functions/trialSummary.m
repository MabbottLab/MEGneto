function trial_summary = trialSummary(paths)

% Descriptiom:
% This function generates a summary of the number of correct trials,
% incorrect trials and total trials per participant.

% How to use:
% call the function with the following command in the command line or as
% presented in the main_template: trialSummary(paths). This function must
% be run until at least after fcp_1_TaskEpoching is run.

% INPUT: paths
%
% OUTPUT: a table structure with the rows representing participants and
% columns representing the total number of trials, total number of correct
% trials, and total number of incorrect trials. 

%% SETUP
% load config JSON with analysis parameters
config      = load_config(paths, paths.name);
config      = config.config;
step        = 'fcp1';
subj_match = ds_pid_match(paths,step);

% initialize output file
trialSummary_output.ppt_trial_summary = 'ppt_trial_summary.mat';

% initialize table to store results
trial_summary = array2table(zeros(height(subj_match), 3));      

%% GENERATING TABLE DATA
for ss = 1:height(subj_match) % for each participant
    dataset = [paths.rawdata '/' subj_match.ds{ss}];
    t0marker = config.task.trialdef.markers.t0marker;
    
    % extract num of total trials
    eventslist = ft_read_event(dataset);
    total_trials = sum(ismember({eventslist.type}, t0marker));
    trial_summary(ss, 1) = {total_trials};
    
    % extract num of correct trials (# 'Correct' in eventslist)
    num_correct = sum(ismember({eventslist.type}, 'Correct'));
    trial_summary(ss, 2) = {num_correct};
    
    % extract num of incorrect trials
    num_incorrect = total_trials - num_correct;
    trial_summary(ss, 3) = {num_incorrect};
    
end
trial_summary.Properties.RowNames = subj_match.pid;
trial_summary.Properties.VariableNames = ["total_trials", ...
                                          "num_correct", ...
                                          "num_incorrect"];
                                      
%% SAVING OUTPUT
fprintf('Saving participant trial summary...\n');
save([paths.anout_grp '/' trialSummary_output.ppt_trial_summary],'trial_summary', '-v7.3')
end
                              