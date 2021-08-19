%% Preamble
% The following script outlines the steps taken to run data analysis on the
% OIRM data set for the faces versus scenes visual content project
% conducted in the summer of 2021.

% Please note that the OIRM Free-Viewing branch of the MEGneto pipeline is
% to be used in tandem with the original MEGneto pipeline (i.e., both are
% needed to conduct the analysis outlined below).

% High-level software/data requirements:
% MATLAB, FieldTrip, MEGneto pipeline, OIRM dataset

%% Table of contents for steps involved in data analysis:
% 1. Obtaining clip markers
% 2. MEG pipeline freeviewing fcp_1-fcp_4 (consists of substeps (a,b)
% 3. freeviewing fcp_4_5_ReinsertingTrials
%.4. MEG pipeline freeviewing fcp_5_freqnalaysis
% 5. freeviewing fcp_5_5_analyzeFreqanalysis

%% 1. Obtaining clip markers

% Please navigate to the Extrapolating_clipTimes folder which contains two
% files (generate_markerTimes.m and read_data.m). The read_data.m file
% contains a function which is used in the generate_markerTimes.m script.To
% run this file you will need the OIRM dataset containing MEG data and .psy
% files (these will indicate movie order for a given run of MEG data). You
% will also need a clipTimes.mat file where for each movie there is an
% array that indicates the frame number where a scene change occured.

% Note: keep in mind for conversion purposes that each frame is equal to 
% 20 samples. 
 
% The generate_markerTimes.m script is used to generate a .mat file
% containing the clip marks (also known as clip samples or clip times), the
% movie order for each run of MEG data, and the clippet marks (or
% samples/times). A pseudocode run-down of how this script operates can be
% found at the top of the script when it is opened.

% Output: a *.mat file containing the clip/clippet times (the file's title
% should be something similiar to clipMarkers_allPpts.mat).

%% 2. MEG pipeline steps
%% a) Preparation for freeviewing fcp1 (filling paths, megne2setup, JSON config)
% Prior to starting with the first step of the pipeline, we must prepare
% some data. The OIRM data is analyzed in TWO separate analysis files - one
% is for all Run 1 data, and the other for all Run 2 data. These two
% separate data analyses are stacked after the fcp4_beamforming step, thus
% all participants/runs are inputted together to the fcp5 step.

% Please refer to the OrderedClips_OIRMRun1_Jun23.m and 
% OrderedClips_OIRMRun2_Jun25.m scripts for reference. 

% Participant files in run 1 include:
% OIRM02_MEG085_20161213_Free.ds
% OIRM03_MEG085_20161219_Free.ds
% OIRM04_MEG085_20161219_Free.ds
% OIRM06_MEG085_20161208_Free.ds
% OIRM07_MEG085_20161216_Free.ds
% OIRM10_MEG092_20161110_FreeViewing.ds
% OIRM12_MEG085_20161209_Free.ds
% OIRM13_MEG085_20161220_Free.ds
% OIRM17_MEG085_20161212_Free.ds
% OIRMP02_MEG092_20161031_Free.ds

% Participant files in run 2 include:
% OIRM02_MEG085_20161213_Free2.ds
% OIRM03_MEG085_20161219_Free2.ds
% OIRM04_MEG085_20161219_Free2.ds
% OIRM06_MEG085_20161208_Free2.ds
% OIRM07_MEG085_20161216_Free2.ds
% OIRM10_MEG092_20161110_FreeViewing2.ds
% OIRM12_MEG085_20161209_Free2.ds
% OIRM13_MEG085_20161220_Free2.ds
% OIRMP02_MEG092_20161031_Free2.ds

% Each analysis should be conducted separately, as mentioned, but both
% analyses for run 1 and run 2 data follow identical steps and will thus
% not be distinguished in this explanation, unless needed.

% Firstly, please follow the standard steps of the MEG pipeline to fill
% out the various paths files. Please ensure the MEG path points to the
% MEGcode folder containg the freeviewing version of all the functions in
% the pipeline! 

% Then, run freeviewing_megne2setup (ensure
% the function name in your main template calls the freeviewing functions
% rather than the normal MEG pipeline. Note that minimal changes exist
% between the standard MEG pipeline functions and the freeviewing
% counterparts. Feel free to explore the two cateogires of functions side 
% by side if you are interested as to where changes occur.

% Move on to the JSON config settings where the fields you will likely
% change are: the email field, the sampling rate field (600hz for OIRM),
% the function field under taskFunc (change it too
% @freeviewingTaskTrialFun), various fields under tasktrialdef (change
% include to "clipChange", tEpoch to [-1,3], Correct to "clipChange",
% t0marker to "clipChange"), the template grid reslution in beamforming to
% 0.65, the atlas file path to "mmp", and in the 11th window (for frequency
% analysis) you will change the first field to 1, the time window of
% interest (toi) to -1:0.033:3, the baseline to [-1, -0.5], and the ROIs to
% match the indices of ROIs in the glasser atlas that you care about
% (please see documentation on OneDrive for the ROIs used in summer 2021).

%% b) Running freeviewing fcp_1-fcp_4
% Next, you will follow standard MEG pipeline steps to complete steps fcp 1
% to fcp 3, but you will ensure that you use the free viewing fcp
% functions. Standard methods are also used to complete fcp4 (
% freeviewing version),  however this step has some noteworthy changes 
% compared to the standard beamforming step such as line 165-180 where 
% fiducials are added. 
% Withing fcp1, the freeviewingTaskTrialFun is called. This function was
% modified to extract information genereated in the clip markers extraction
% process mentioned earlier. Essentially, this step ensures that trials are
% defined by clip and clippet changes and re-orders all movie clips to
% follow a consistent format (1-5 for run 1 data and 6-10 for run 2 data).
% Also note that any trials that get rejected are noted in this step so
% they can be re-inserted later on as NaNs. This is important since trials
% are defined by clips/clippets and we want to average data across
% consistent clips/clippets later to ensure the same visual content is
% present across trials.

%% 3. freeviewing fcp_4_5_ReinsertingTrials 
% After fcp4 beamforming is completed we want to stack the run 1 and run 2
% data for each participant that has two runs. We also want to insert blank
% trials on a per-particiapnt basis for trials that were removed in fcp1.
% The freeviewing fcp_4_5 step facilitates this and is separated into
% sections which should be run in order.The script extracts a 
% rejected_trls.mat file located in each participant's folder (there will 
% be one for each run), and insert NaN columns at those specified indices.

% From here on out you should choose where you will be storing the
% remainder of your data (i.e., which analysis folder and thus which paths
% variable will be loaded, as your runs of data are now combined, 
% not separated)

%% 4. freeviewing fcp_5_freqanalysis
% Now we will run the frequency analysis step. The configuration for this
% step was completed in the JSON file set up so we do not need to worry
% about that. However we should pay attention to the following lines 
% of code:
%%% line 70 - 71: load the feature vector that designates faces versus other
% visual content (faces = 1, other = 0).
%%% lines 103 - 114: identify NaN trials (the ones we reinserted in the step
% above) and remove them since this messes up stats analyses later on. Note
% that because the analysis in the end doesnt separate out individual
% trials, and instead AVERAGES across all trials containing a specific
% category of visual stimulus, we did not need to reinsert the NaN trials
% but it is good that we have this option in case we want to do analyses on
% a trial by trial basis. 
%%% line 121-123:  this is where we specify what trials we are including. 
% To use all trials (and thus not separate by category of visual content)
% leave line 123 uncommented, else you must comment line 123 and uncomment
% line 122 where you specify 1 for faces and 0 for other content. 
%%% lines 147 - 156:this is where we save the output. Be sure that you save
% your output in a location that makes sense and under a name that is
% sensisble (e.g. include "faces" in the name if you are including only
% trials with faces). 

% The output of this function is a power spectrum with the following
% dimensions: participants x ROIs x frequency x time. 

%% 5. freeviewing fcp_5_5_analyzeFreqanalysis 
% Finally we are at the point where we will run a statistical analysis to
% test hypotheses and generate figures! 

% This script is separated into sections which you should run in order. 
% A brief summary is as follows:
%%% The first section of the analysis is for the first aim: Is there a
% change in neural oscillatory power after stimulus presentation? The data
% for this is not separated by cateogory of visual stimulus. First we
% prepare data for stats input (get a grand avg across all ppts, get a
% grand avg pre-stimulus, get a grand avg post-stimulus). Then we perform a
% monte carlo permutation to shuffle groups labels (pre stimulus, post
% stimulus) and create a null distribution of t-statistics. Then we use a
% one way t-test to determine whether there is a significant difference 
% post stimulus comapred to pre stimulus. Following this, we generate
% several plots to communicate information about the hypotheses.

%%% The second section of the analysis is for the second aim: Is there a
% change in neural oscillatory power for visual stimuli containing faces
% compared to scenes? First we must prepare the data by averaging over
% aprticiapnt for the faces category and then for the scenes catoegyr.
% Again we do a monte carlo permutation followed by a one way t-test to
% determine if there is an increase in power for faces compared to scenes.
% Then, we generate severall figures to represent this information. Please
% note that a spreadhseet in OneDrive in the analysis notes folder 
% called FacesVScenes_figure can be found which summarizes these findings 
% as well. 





