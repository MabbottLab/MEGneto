%% NBS FOR LOOP run per frequency band and contrast for t-threshold values of 0.05 uptil 4.6

%STEP 1: Start NBS

%STEP 2: Select design matrix, frequency band datamatrix and contrats of choice

%STEP 3: Set the t-threshold to the lowest value of 0.05

%STEP 4: Run NBS

%STEP 5: Save NBS result under any name using the drop down menu in the GUI 

%STEP 6: PUll up NBS result file into workspace and open "nbs" variable
%within the file in the workspace

%STEP 7: Run this loop
%% NBS LOOP
thresh_values = [0.05:0.01:0.1 0.2:0.1:4.6]; %set threshold values to test
%T-thresholds from 0.05 will increase by increments of 0.01 until 0.1
%From 0.1 thresholds will increase by increments of 0.1 until 4.6
p_vals = []; %Where p-values for each iteration will be saved if greater than 0.05

for thr = 1:length(thresh_values)
    nbs.UI.thresh.ui = num2str(thresh_values(thr));
    NBSrun(nbs.UI,[])
    global nbs
    if isempty(nbs.NBS.pval)
        p_vals(thr) = NaN;%Save thresholds that are not significant as NAN in p_vals variable
    else
        p_vals(thr) = nbs.NBS.pval;
    end
end

%
%% IF NO NBS POP-UP WINDOWS APPEAR, NO SIGNIFICANT CONNECTIONS WERE FOUND
%PROCEED TO RUN NBS AND REPEAT STEPS 1-7 WITH THE NEXT CONTRAST OR
%FREQUENCY BAND OF CHOICE
%save NBS Results from step 5 under same name, it will replace the one in
%the workspace.
%% IF NBS POP-UP WINDOWS APPEAR
% STEP 9: close all NBS results tabs 

% STEP 10: find out which threshold gave the lowest p_vals
[P,I] = min (p_vals);
Low_T = thresh_values(I);
Low_T
% STEP 9: Run NBS in command line again to pull up GUI

% STEP 11: Insert the lowest threshold value (Low_T) printed on command
% line or check variable in workspace

% STEP 12: Run NBS with that threshold and save file as
% AnalysisName_NBSResults_freq.mat

%Repeat for each contrast for each frequency band*****