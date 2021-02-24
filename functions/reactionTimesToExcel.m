function reactionTimesToExcel(reactionTimeSummary, outdir)

%% Purpose:
% This function saves the reaction time summary data (output of the
% reactionTimes function) as an excel file for each participant. This
% function is to be run only after the reactionTimes function has been run
% and is only relevant if the user wishes to generate excel files
% containing this data.

%% Inputs:
% 1. reactionTimeSummary
% variable containing the output of the reactionTimes
% function. If you click "open" then go into the "group" folder within the
% "analysis" folder of your analysis (specified by what you set 
% "analysis_name" to be on the main template), and double click the
% 'reactionTimeSummary.mat' file, you will see the variable
% reactionTimeSummary. If you follow those steps then the variable you pass
% in for this input will indeed be reactionTimeSummary.

% 2. outdir
% a string specifying which folder the outputs of this function should be
% stored in (e.g. '/home/your_username' for your home folder)

%% Outputs:
% Excel files for each participant containing their reaction time summary
% data (as extracted from the output of the reactionTimes function)

%% Code

for i = 1:height(reactionTimeSummary)
    filename = strcat('reactionTimeSummary_', reactionTimeSummary.Properties.RowNames{i}, '.xlsx');
    fullpath = fullfile(outdir, filename);
    writetable(reactionTimeSummary.data{i,1}{1,1}, fullpath)
end 

end 