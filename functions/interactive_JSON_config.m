function interactive_JSON_config(paths, megneto_path)

% This user-interactive file aids users in filling out the JSON config file
% which specifies parameters and instructions for various parts of the
% pipeline. This step is done just before fcp_1 is run. 

% For more help and documentation on the JSON config file, the user can
% consult the config file documentation on the Mabbott Lab GitHub.

% Inputs: 
%       - path to empty template JSON config file
%       - paths to analysis and project folder

% Outputs: populated JSON config file

% Notes:
%       -str2double is used to convert string inputs to numbers, where needed
%       -str2num is used to convert string inputs to a Nx1 double where N > 1
%       -split is used to separate strings into cells 
%       -strcmp is used to convert a string input into a logical input

% The overall flow of this code is to go section by section of the template
% JSON config file and set up an input dialog box with pre-entered values
% that are default values/show the user the proper input format. After each
% section, that portion of the JSON config file is filled in. After all
% sections are completed, the JSOn config file is encoded to its JSON
% format and formatting (spacing, tabbing, etc.) is conducted.

%% Load in the blank config template
empty_config             = fullfile(megneto_path, '/configs/empty_config.json');
decode                   = loadjson(empty_config);

%% Setting up information for the dialog boxes

% prompts for the dialog boxes
prompt{1}     = {'Enter contact email', ...
                 'Enter epoching period for resting state data',... 
                 'Enter threshold for headmotion'};
prompt{2}     = {'Enter 0 (no) or 1 (yes) to indicate if artifact detection is desired',...
                 'Enter yes or no if you wish to bandpass filter for muscle artifact '...
                 'Enter the bandpass frequency for muscle artifact [x y]',...
                 'Enter the filter order for the muscle artifact cleaning',...
                 'Enter the filter type the muscle artifact (e.g. but for butterworth)'...
                 'Enter yes or no to indicate if a hilbert transform should be done'...
                 'Enter a window length (in terms of time) for the boxcar/moving average filter'...
                 'Enter the cut off frequency for muscle artifact cleaning'...
                 'Enter desired time length for trial padding'};
prompt{3}     = {'Enter desired time length for filter padding'...
                 'Enter desired time length for artifact padding'...
                 'Enter cutoff frequency for jump artifact detection'...
                 'Enter 0 (no) or 1 (yes) to indicate if ICA should be done'...
                 'Enter 0 (no) or 1 (yes) to indicate if noisy trials should be removed'...
                 'Enter 0 (no) or 1 (yes) to indicate if bad trials should be removed'};
prompt{4}     = {'Specify which data channel to analyse for MEG data',...
                 'Enter yes or no to indicate if a notch filter for line noise should be applied'...
                 'Enter the frequency of the line noise [x, y]',...
                 'Enter yes or no to indicate if a bandpass filter should be applied'...
                 'Enter bandpass frquency [x, y]',...
                 'Enter the filter order', ...
                 'Enter the sampling rate of the data',...
                 'Enter the MEG model'};
prompt{5}     = {'Enter the name of the function used to parse data into trials',...
                 'Specify the type of function for the previous filed',...
                 'Enter file',... 
                 'Enter within_file_path'}; 
prompt{6}     = {'Enter 0 (no) or 1 (yes) to specify if there is rest data', ...
                 'Enter the name of marker to epoch around',...
                 'Enter marker to gain reaction time information on', ... % need clarity
                 'Enter  true or false to specify whether you only want information on the number of trials'...
                 'Specify delay time to correct for',... % need clarity
                 'Enter epoching period',...
                 'Specify markers that distinguish a correct trial',...
                 'Specify the presentation stimulus marker'};
prompt{7}     = {'Specify the headmodel form',...
                 'Specify the headmodel units',...
                 'Specify the template MNI grid resolution in cm',...
                 'Specify yes or no for tight grid',... % need clarity
                 'Specify how much the innermost surface should be moved inward to constrain sources to be inside the source compartment',...
                 'Specify the coordinate system',... 
                 'Specify the atlas filepath',...
                 'Specify the input coordinate',... % need clarity
                 'Specify the plotting method of MRI volumes',...
                 'Specify the dimension to slice on if the plotting field method is "slice"'};
prompt{8}     = {'Specify the number of slices if the plotting field method is "slice"',...
                 'Yes or no to specify whether we want to warp the model to the T1 model',...
                 'Enter yes or no to indicate whether non-linear normalization should be used.',...
                 'Specify grid units',...
                 'Enter yes or no to specify whether we want to normalize, which addresses depth bias',...
                 'Enter yes or no if a covariance matrix should be computed',...
                 'Specify the window length for covariance matrix computation',...
                 'Enter the vartrllength',... % need clarity
                 'Enter yes or no to specify if data should be projected to dominant eigenvector'};
prompt{9}     = {'Enter yes or no to specify if trials should be separated',...
                 'Enter yes or no if the filter should be kept for projection'...
                 'Enter yes or no to construct a filter from single trials and apply it to single trials',...
                 'Enter beamforming method',...
                 'Enter method of generating a representative timeseries'};
prompt{10}    = {'Enter the metric for connectivity analysis',...
                 'Specify the frequency bands',...
                 'Specify the frequency band  names, in order of the previous field',...
                 'Specify the method of generating a representative timeseries for each ROI'};

% default inputs for the dialog boxes
definput{1}   = {'firstname.lastname@sickkids.ca', '30', '10'};
definput{2}   = {'1', 'yes','[110,140]', '8', 'but', 'yes', '0.2', '30',...
                '0.5'};
definput{3}   = {'0.1', '0.1', '35', '1', '1', '1'};
definput{4}   = {'MEG,MEGREF,REFGRAD,REFMAG','yes', '[60,120]', 'yes',...
                 '[1,150]', '5', '300', 'CTF151.lay'};
definput{5}   = {'@searchTaskTrialFun', 'anonymous', '', '__base_function'};
definput{6}   = {'0', 'Correct', 'Correct', 'false', '0.023',...
                 '[-2.0, 2.0]', 'LeftCorrect,RightCorrect',...
                 'OfflineLightOn'};
definput{7}   = {'singleshell', 'cm', '1', 'yes', '-0.8', 'spm',...
                 '/template/atlas/aal/ROI_MNI_V4.nii', 'mni', 'slice', '2'};
definput{8}   = {'20', 'yes', 'yes', 'cm', 'no', 'yes', 'all', '2', 'yes'};
definput{9}   = {'yes', 'yes', 'yes', 'lcmv', 'mean'};
definput{10}  = {'wpli_debiased', '[4,7;8,12;13,29;30,59;60,100]',...
                 'theta,alpha,beta,lowgamma,highgamma', 'max'};

% titles of dialog boxes
dlg_title{1}  = 'Part 1 of 10: config.contact and config.epoching';
dlg_title{2}  = 'Part 2 of 10: config.cleaningoptions';
dlg_title{3}  = 'Part 3 of 10: config.cleaningoptions';
dlg_title{4}  = 'Part 4 of 10: config.filteringParameters';
dlg_title{5}  = 'Part 5 of 10: config.taskFunc';
dlg_title{6}  = 'Part 6 of 10: config.task';
dlg_title{7}  = 'Part 7 of 10: config.beamforming';
dlg_title{8}  = 'Part 8 of 10: config.beamforming';
dlg_title{9}  = 'Part 9 of 10: config.beamforming';
dlg_title{10} = 'Part 10 of 10: config.connectivity';

% dimensions for dialog box
dims          = [1, 150];

%% Generating dialog boxes (10 total)

for cfg_part = 1:10
    answers = inputdlg(prompt{cfg_part}, dlg_title{cfg_part}, dims, definput{cfg_part});
    switch cfg_part
        case 1
            decode.config.contact = {answers{1}};
            decode.config.epoching.period = str2double(answers{2}); 
            decode.config.epoching.headMotion.thr = str2double(answers{3});
        case 2
            decode.config.cleaningOptions.artifact.detection = str2double(answers{1});
            decode.config.cleaningOptions.artifact.muscle.bpfilter = answers{2};
            decode.config.cleaningOptions.artifact.muscle.bpfreq = str2num(answers{3});
            decode.config.cleaningOptions.artifact.muscle.bpfiltord = str2double(answers{4});
            decode.config.cleaningOptions.artifact.muscle.bpfilttype = answers{5};
            decode.config.cleaningOptions.artifact.muscle.hilbert = answers{6};
            decode.config.cleaningOptions.artifact.muscle.boxcar = str2double(answers{7});
            decode.config.cleaningOptions.artifact.muscle.cutoff = str2double(answers{8});
            decode.config.cleaningOptions.artifact.muscle.trlpadding = str2double(answers{9});
        case 3
            decode.config.cleaningOptions.artifact.muscle.fltpadding = str2double(answers{1});
            decode.config.cleaningOptions.artifact.muscle.artpadding = str2double(answers{2});
            decode.config.cleaningOptions.artifact.jump.cutoff = str2double(answers{3});
            decode.config.cleaningOptions.artifact.icaClean = str2double(answers{4});
            decode.config.cleaningOptions.artifact.rmNoisyTrials = str2double(answers{5});
            decode.config.cleaningOptions.rmBadChannels = str2double(answers{6});
        case 4
            decode.config.filteringParameters.channel = split(answers{1}, ',');
            decode.config.filteringParameters.dftfilter = answers{2};
            decode.config.filteringParameters.dftfreq = str2num(answers{3});
            decode.config.filteringParameters.bpfilter = answers{4};
            decode.config.filteringParameters.bpfreq = str2num(answers{5});
            decode.config.filteringParameters.bpfiltord = str2double(answers{6});
            decode.config.filteringParameters.sampleRate = str2double(answers{7});
            decode.config.filteringParameters.CTFlayout = answers{8};
        case 5
            decode.config.taskFunc.function = answers{1};
            decode.config.taskFunc.type = answers{2};
            decode.config.taskFunc.file = answers{3};
            decode.config.taskFunc.within_file_path = answers{4};
        case 6
            decode.config.task.isRest = str2double(answers{1});
            decode.config.task.trialdef.details.name = {answers{2}};
            decode.config.task.trialdef.details.include = {answers{3}};
            decode.config.task.trialdef.details.countOnly = strcmp(answers{4}, 'true');
            decode.config.task.trialdef.parameters.t0shift = str2double(answers{5});
            decode.config.task.trialdef.parameters.tEpoch = str2num(answers{6});
            decode.config.task.trialdef.markers.Correct = split(answers{7}, ',');
            decode.config.task.trialdef.markers.t0marker = answers{8};
        case 7
            decode.config.beamforming.headmodel.method = answers{1};
            decode.config.beamforming.headmodel.units = answers{2};
            decode.config.beamforming.template.grid.resolution = str2double(answers{3});
            decode.config.beamforming.template.grid.tight = answers{4};
            decode.config.beamforming.template.grid.inwardshift = str2double(answers{5});
            decode.config.beamforming.template.coordsys = answers{6};
            decode.config.beamforming.atlas.filepath = answers{7};
            decode.config.beamforming.atlas.inputcoord = answers{8};
            decode.config.beamforming.checkMRIvolumes.method = answers{9};
            decode.config.beamforming.checkMRIvolumes.slicesdim = str2double(answers{10});
        case 8
            decode.config.beamforming.checkMRIvolumes.nslices = str2double(answers{1});
            decode.config.beamforming.subj.grid.warpmni = answers{2};
            decode.config.beamforming.subj.grid.nonlinear = answers{3};
            decode.config.beamforming.subj.grid.unit = answers{4};
            decode.config.beamforming.leadfield.normalize = answers{5};
            decode.config.beamforming.timeDomain.covariance = answers{6};
            decode.config.beamforming.timeDomain.covariancewindow = answers{7};
            decode.config.beamforming.timeDomain.vartrllength = str2double(answers{8});
            decode.config.beamforming.timeDomain.projectmom = answers{9};
        case 9
            decode.config.beamforming.options.keeptrials = answers{1};
            decode.config.beamforming.options.keepfilter = answers{2};
            decode.config.beamforming.options.rawtrial = answers{3};
            decode.config.beamforming.method = answers{4};
            decode.config.beamforming.rep_timeseries = answers{5};
        case 10
            decode.config.connectivity.method = answers{1};
            decode.config.connectivity.filt_freqs = str2num(answers{2});
            decode.config.connectivity.freq_names = split(answers{3}, ',');
            decode.config.connectivity.collapse_band = answers{4};
    end 
            
end

%% Save JSON config file

save_to_json(decode, paths.mainconf, true);