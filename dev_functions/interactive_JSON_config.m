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
empty_config = fullfile(megneto_path, '/configs/empty_config.json');
jsonText = fileread(empty_config);
decode = jsondecode(jsonText);

%% Filling in sections of the config file
%% config.contact and config.epoching
contact_epoching_prompt = {'Enter contact email', ...
'Enter epoching period for resting state data',... 
'Enter threshold for headmotion'};

contact_epoching_dlgtitle = 'config.contact and config.epoching';
dims = [1, 70];

contact_epoching_definput = {'firstname.lastname@sickkids.ca', '30', '10'};
        
contact_epoching_answer = inputdlg(contact_epoching_prompt, contact_epoching_dlgtitle, dims, contact_epoching_definput);

decode.config.contact = {contact_epoching_answer{1}};
decode.config.epoching.period = str2double(contact_epoching_answer{2}); 
decode.config.epoching.headMotion.thr = str2double(contact_epoching_answer{3});

%% config.cleaningoptions 
cleaningoptions_prompt = {'Enter 0 (no) or 1 (yes) to indicate if artifact detection is desired',...
'Enter yes or no if you wish to bandpass filter for muscle artifact '...
'Enter the bandpass frequency for muscle artifact [x y]',...
'Enter the filter order for the muscle artifact cleaning',...
'Enter the filter type the muscle artifact (e.g. but for butterworth)'...
'Enter yes or no to indicate if a hilbert transform should be done'...
'Enter a window length (in terms of time) for the boxcar/moving average filter'...
'Enter the cut off frequency for muscle artifact cleaning'...
'Enter desired time length for trial padding'};

cleaningoptions_dlgtitle = 'config.cleaningoptions';
dims = [1, 70];

cleaningoptions_definput = {'1', 'yes','[110;140]', '8', 'but', 'yes', '0.2',...
            '30', '0.5'};
        
cleaningoptions_answer = inputdlg(cleaningoptions_prompt, cleaningoptions_dlgtitle, dims, cleaningoptions_definput);

decode.config.cleaningOptions.artifact.detection = str2double(cleaningoptions_answer{1});
decode.config.cleaningOptions.artifact.muscle.bpfilter = cleaningoptions_answer{2};
decode.config.cleaningOptions.artifact.muscle.bpfreq = str2num(cleaningoptions_answer{3});
decode.config.cleaningOptions.artifact.muscle.bpfiltord = str2double(cleaningoptions_answer{4});
decode.config.cleaningOptions.artifact.muscle.bpfilttype = cleaningoptions_answer{5};
decode.config.cleaningOptions.artifact.muscle.hilbert = cleaningoptions_answer{6};
decode.config.cleaningOptions.artifact.muscle.boxcar = str2double(cleaningoptions_answer{7});
decode.config.cleaningOptions.artifact.muscle.cutoff = str2double(cleaningoptions_answer{8});
decode.config.cleaningOptions.artifact.muscle.trlpadding = str2double(cleaningoptions_answer{9});

%% config.cleaningoptions part 2
cleaningoptions2_prompt = {'Enter desired time length for filter padding'...
'Enter desired time length for artifact padding'...
'Enter cutoff frequency for jump artifact detection'...
'Enter 0 (no) or 1 (yes) to indicate if ICA should be done'...
'Enter 0 (no) or 1 (yes) to indicate if noisy trials should be removed'...
'Enter 0 (no) or 1 (yes) to indicate if bad trials should be removed'};

cleaningoptions2_dlgtitle = 'config.cleaningoptions';
dims = [1, 70];

cleaningoptions2_definput = {'0.1', '0.1', '35', '1', '1', '1'};
        
cleaningoptions2_answer = inputdlg(cleaningoptions2_prompt, cleaningoptions2_dlgtitle, dims, cleaningoptions2_definput);

decode.config.cleaningOptions.artifact.muscle.fltpadding = str2double(cleaningoptions2_answer{1});
decode.config.cleaningOptions.artifact.muscle.artpadding = str2double(cleaningoptions2_answer{2});
decode.config.cleaningOptions.artifact.jump.cutoff = str2double(cleaningoptions2_answer{3});
decode.config.cleaningOptions.artifact.icaClean = str2double(cleaningoptions2_answer{4});
decode.config.cleaningOptions.artifact.rmNoisyTrials = str2double(cleaningoptions2_answer{5});
decode.config.cleaningOptions.rmBadChannels = str2double(cleaningoptions2_answer{6});

%% config.filteringParameters
filteringParams_prompt = {'Specify which data channel to analyse for MEG data',...
 'Enter yes or no to indicate if a notch filter for line noise should be applied'...
 'Enter thw frequency of the line noise [x, y]',...
 'Enter yes or no to indicate if a bandpass filter should be applied'...
 'Enter bandpass frquency [x, y]',...
 'Enter the filter order', ...
 'Enter the sampling rate of the data',...
 'Enter the MEG model'};

filteringParams_dlgtitle = 'config.filteringParameters';
dims = [1, 70];

filteringParams_definput = {'MEG, MEGREF, REFGRAD, REFMAG',...
'yes', '[60;120]', 'yes', '[1;150]', '5', '300', 'CTF151.lay'};

filteringParams_answer = inputdlg(filteringParams_prompt, filteringParams_dlgtitle, dims, filteringParams_definput);

decode.config.filteringParameters.channel = split(filteringParams_answer{1}, ',');
decode.config.filteringParameters.dftfilter = filteringParams_answer{2};
decode.config.filteringParameters.dftfreq = str2num(filteringParams_answer{3});
decode.config.filteringParameters.bpfilter = filteringParams_answer{4};
decode.config.filteringParameters.bpfreq = str2num(filteringParams_answer{5});
decode.config.filteringParameters.bpfiltord = str2double(filteringParams_answer{6});
decode.config.filteringParameters.sampleRate = str2double(filteringParams_answer{7});
decode.config.filteringParameters.CTFlayout = filteringParams_answer{8};

%% config.taskFunc
taskFunc_prompt = {'Enter the name of the function used to parse data into trials',...
'Specify the type of function for the previous filed',...
'Enter file',... % ISSUE: unclear what this field is
'Enter within_file_path'}; % ISSUE: unclear what this field is 
% ISSUE: removed "workspace" field as it is unclear what it is & due to its
% struct nature

taskFunc_dlgtitle = 'config.taskFunc';
dims = [1, 70];

taskFunc_definput = {'@searchTaskTrialFun', 'anonymous', '', '__base_function'};
        
taskFunc_answer = inputdlg(taskFunc_prompt, taskFunc_dlgtitle, dims, taskFunc_definput); 

decode.config.taskFunc.function = taskFunc_answer{1};
decode.config.taskFunc.type = taskFunc_answer{2};
decode.config.taskFunc.file = taskFunc_answer{3};
decode.config.taskFunc.within_file_path = taskFunc_answer{4};

%% config.task
task_prompt = {'Enter 0 (no) or 1 (yes) to specify if there is rest data', ...
'Enter the name of marker to epoch around',...
'Enter the name of marker to include once (leave as blank)',... % need clarity
'Enter name of marker to exclude (leave as blank)',... % need clarity
'Enter marker to gain reaction time information on', ... % need clarity
'Enter  true or false to specify whether you only want information on the number of trials'...
'Specify delay time to correct for',... % need clarity
'Enter epoching period',...
'Specify markers that distinguish a correct trial',...
'Specify the presentation stimulus marker'};

task_dlgtitle = 'config.task';
dims = [1, 70];

task_definput = {'0', 'Correct', '', '', 'Correct', 'false',...
'0.023', '[-2.0;2.0]', 'LeftCorrect, RightCorrect', 'OfflineLightOn'};
        
task_answer = inputdlg(task_prompt, task_dlgtitle, dims, task_definput);

decode.config.task.isRest = str2double(task_answer{1});
decode.config.task.trialdef.details.name = {task_answer{2}};
decode.config.task.trialdef.details.includeOnce = {task_answer{3}};
decode.config.task.trialdef.details.exclude = {task_answer{4}};
decode.config.task.trialdef.details.include = {task_answer{5}};
decode.config.task.trialdef.details.countOnly = strcmp(task_answer{6}, 'true');
decode.config.task.trialdef.parameters.t0shift = str2double(task_answer{7});
decode.config.task.trialdef.parameters.tEpoch = str2num(task_answer{8});
decode.config.task.trialdef.markers.Correct = split(task_answer{9}, ',');
decode.config.task.trialdef.markers.t0marker = task_answer{10};

%% config.beamformin[firstname.lastname@sickkids.ca]g
beamforming_prompt = {'Specify the headmodel form',...
'Specify the headmodel units',...
'Specify the template MNI grid resolution in mm',...
'Specify yes or no for tight grid',... % need clarity
'Specify how much the innermost surface should be moved inward to constrain sources to be inside the source compartment',...
'Specify the coordinate system',... 
'Specify the atlas filepath',...
'Specify the input coordinate',... % need clarity
'Specify the plotting method of MRI volumes',...
'Specify the dimension to slice on if the plotting field method is "slice"'};

beamforming_dlgtitle = 'config.beamforming';
dims = [1, 70];

beamforming_definput = {'singleshell', 'cm', '1', 'yes', '-0.8', 'spm',... 
'/template/atlas/aal/ROI_MNI_V4.nii', 'mni', 'slice', '2'};
        
beamforming_answer = inputdlg(beamforming_prompt, beamforming_dlgtitle, dims, beamforming_definput); 

decode.config.beamforming.headmodel.method = beamforming_answer{1};
decode.config.beamforming.headmodel.units = beamforming_answer{2};
decode.config.beamforming.template.grid.resolution = str2double(beamforming_answer{3});
decode.config.beamforming.template.grid.tight = beamforming_answer{4};
decode.config.beamforming.template.grid.inwardshift = str2double(beamforming_answer{5});
decode.config.beamforming.template.coordsys = beamforming_answer{6};
decode.config.beamforming.atlas.filepath = beamforming_answer{7};
decode.config.beamforming.atlas.inputcoord = beamforming_answer{8};
decode.config.beamforming.checkMRIvolumes.method = beamforming_answer{9};
decode.config.beamforming.checkMRIvolumes.slicesdim = str2double(beamforming_answer{10});

%% config.beamforming part 2
beamforming2_prompt = {'Specify the number of slices if the plotting field method is "slice"',...
'Yes or no to specify whether we want to warp the model to the T1 model',...
'Enter yes or no to indicate whether non-linear normalization should be used.',...
'Specify grid units',...
'Enter yes or no to specify whether we want to normalize, which addresses depth bias',...
'Enter yes or no if a covariance matrix should be computed',...
'Specify the window length for covariance matrix computation',...
'Enter vartrllength',... % need to remove from other templates and here!!
'Enter yes or no to specify if data should be projected to dominant eigenvector'};

beamforming2_dlgtitle = 'config.beamforming';
dims = [1, 70];

beamforming2_definput = {'20', 'yes', 'yes', 'cm', 'no', 'yes', 'all', '2', 'yes'};
        
beamforming2_answer = inputdlg(beamforming2_prompt, beamforming2_dlgtitle, dims, beamforming2_definput); 

decode.config.beamforming.checkMRIvolumes.nslices = str2double(beamforming2_answer{1});
decode.config.beamforming.subj.grid.warpmni = beamforming2_answer{2};
decode.config.beamforming.subj.grid.nonlinear = beamforming2_answer{3};
decode.config.beamforming.subj.grid.unit = beamforming2_answer{4};
decode.config.beamforming.leadfield.normalize = beamforming2_answer{5};
decode.config.beamforming.timeDomain.covariance = beamforming2_answer{6};
decode.config.beamforming.timeDomain.covariancewindow = beamforming2_answer{7};
decode.config.beamforming.timeDomain.vartrllength = str2double(beamforming2_answer{8});
decode.config.beamforming.timeDomain.projectmom = beamforming2_answer{9};

%% config.beamforming part 3
beamforming3_prompt = {'Enter yes or no to specify if trials should be separated',...
'Enter yes or no if the filter should be kept for projection'...
'Enter yes or no',...
'Enter beamforming method',...
'Enter method of generating a representative timeseries'};

beamforming3_dlgtitle = 'config.beamforming';
dims = [1, 70];

beamforming3_definput = {'yes', 'yes', 'yes', 'lcmv', 'mean'};
        
beamforming3_answer = inputdlg(beamforming3_prompt, beamforming3_dlgtitle, dims, beamforming3_definput); 

decode.config.beamforming.options.keeptrials = beamforming3_answer{1};
decode.config.beamforming.options.keepfilter = beamforming3_answer{2};
decode.config.beamforming.options.rawtrial = beamforming3_answer{3};
decode.config.beamforming.method = beamforming3_answer{4};
decode.config.beamforming.rep_timeseries = beamforming3_answer{5};

%% config.connectivity
connectivity_prompt = {'Enter the metric for connectivity analysis',...
'Specify the frequency bands',...
'Specify the method of generating a representative timeseries for each ROI'};

connectivity_dlgtitle = 'config.connectivity';
dims = [1, 70];

connectivity_definput = {'wpli_debiased', '[4,7;8,12;13,29;30,59;60,100]', 'max'};
connectivity_answer = inputdlg(connectivity_prompt, connectivity_dlgtitle, dims, connectivity_definput); 

decode.config.connectivity.method = connectivity_answer{1};
decode.config.connectivity.filt_freqs = str2num(connectivity_answer{2});
decode.config.connectivity.collapse_band = connectivity_answer{3};

%% encode the json into proper format 
encode = jsonencode(decode);

%% add spacing and indendation 
encode = strrep(encode, '{}', sprintf('[{}]')); % for taskFunc.workspace (untouched by user) 
encode = strrep(encode, '[{', sprintf('[\n\t\t\t\t{'));
encode = strrep(encode, '}]', sprintf('}\n\t\t\t]'));
encode = strrep(encode, '[[', sprintf('[\n\t['));
encode = strrep(encode, ']]', sprintf(']\n\t]'));
encode = strrep(encode, '"," ', sprintf('", "'));
encode = strrep(encode, ':{', sprintf(':{\n\t'));
encode = strrep(encode, '],', sprintf('],\n\t'));
encode = strrep(encode, '},', sprintf('},\n\t'));
encode = strrep(encode, ':', '  : ');
encode = strrep(encode, ',"', sprintf(',\n\t"'));
encode = strrep(encode, '{"', sprintf('{\n\t"'));

%% tabbing specific fields
%% section 1: contact/epoching
encode = strrep(encode, '"contact', sprintf('\t"contact'));
encode = strrep(encode, '"epoching', sprintf('\t"epoching'));
encode = strrep(encode, '"period', sprintf('\t\t"period'));
encode = strrep(encode, '"headMotion', sprintf('\t\t"headMotion'));
encode = strrep(encode, '"thr', sprintf('\t\t\t"thr'));

%% section 2: cleaningOptions
encode = strrep(encode, '"cleaningOptions', sprintf('\t"cleaningOptions'));
encode = strrep(encode, '"artifact', sprintf('\t\t\"artifact'));
encode = strrep(encode, '"detection', sprintf('\t\t\t"detection'));
encode = strrep(encode, '"muscle', sprintf('\t\t\t"muscle'));
encode = strrep(encode, '"bp', sprintf('\t\t\t\t"bp'));
encode = strrep(encode, '"hilber', sprintf('\t\t\t\t"hilbert'));
encode = strrep(encode, '"boxcar', sprintf('\t\t\t\t"boxcar'));
encode = strrep(encode, '"cutoff', sprintf('\t\t\t\t"cutoff'));
encode = strrep(encode, '"trl', sprintf('\t\t\t\t"trl'));
encode = strrep(encode, '"flt', sprintf('\t\t\t\t"flt'));
encode = strrep(encode, '"artpadding', sprintf('\t\t\t\t"artpadding'));
encode = strrep(encode, '"jump', sprintf('\t\t\t"jump'));
encode = strrep(encode, '"icaClean', sprintf('\t\t\t"icaClean'));
%fid = fopen('Testing.json', 'w');
encode = strrep(encode, '"rmNoisyTrials', sprintf('\t\t\t"rmNoisyTrials'));
encode = strrep(encode, '"rmBadChannels', sprintf('\t\t"rmBadChannels'));

%fid = fopen('Testing.json', 'w');
%% section 3: filteringParameters
encode = strrep(encode, '"filteringParameters', sprintf('\t"filteringParameters'));
encode = strrep(encode, '"channel', sprintf('\t\t\t\t"channel'));
%fid = fopen('Testing.json', 'w');
encode = strrep(encode, '"dft', sprintf('\t\t\t\t"dft'));
encode = strrep(encode, '"sampleRate', sprintf('\t\t\t\t"sampleRate'));
encode = strrep(encode, '"CTFlayout', sprintf('\t\t\t\t"CTFlayout'));
%fid = fopen('Testing.json', 'w'););

%% section 4: taskFunc
encode = strrep(encode, '"taskFunc', sprintf('\t"taskFunc'));
encode = strrep(encode, '"function', sprintf('\t\t"function'));
encode = strrep(encode, '"type', sprintf('\t\t"type'));
encode = strrep(encode, '"file"', sprintf('\t\t"file"'));
encode = strrep(encode, '"workspace', sprintf('\t\t"workspace'));
encode = strrep(encode, '"within_file_path', sprintf('\t\t"within_file_path'));

%% section 5: task
encode = strrep(encode, '"task"', sprintf('\t"task"'));
encode = strrep(encode, '"isRest', sprintf('\t\t"isRest'));
encode = strrep(encode, '"trialdef', sprintf('\t\t"trialdef'));
encode = strrep(encode, '"details', sprintf('\t\t\t"details'));
encode = strrep(encode, '"name', sprintf('\t\t\t\t"name'));
encode = strrep(encode, '"include', sprintf('\t\t\t\t"include'));
encode = strrep(encode, '"exclude', sprintf('\t\t\t\t"exclude'));
encode = strrep(encode, '"count', sprintf('\t\t\t\t"count'));
encode = strrep(encode, '"parameters', sprintf('\t\t\t"parameters'));
encode = strrep(encode, '"t0shift', sprintf('\t\t\t\t"t0shift'));
encode = strrep(encode, '"tEpoch', sprintf('\t\t\t\t"tEpoch'));
encode = strrep(encode, '"markers', sprintf('\t\t\t"markers'));
encode = strrep(encode, '"Correct"  :', sprintf('\t\t\t\t"Correct"  :'));
encode = strrep(encode, '"t0marker"', sprintf('\t\t\t\t"t0marker"'));

%% section 6: beamforming
encode = strrep(encode, '"beamforming', sprintf('\t"beamforming'));
encode = strrep(encode, '"headmodel', sprintf('\t\t"headmodel'));
encode = strrep(encode, '"method', sprintf('\t\t"method')); % ISSUE %
encode = strrep(encode, '"units"', sprintf('\t\t\t"units"'));
encode = strrep(encode, '"template', sprintf('\t\t"template'));
encode = strrep(encode, '"grid', sprintf('\t\t\t"grid'));
encode = strrep(encode, '"resolution', sprintf('\t\t\t\t"resolution'));
encode = strrep(encode, '"tight', sprintf('\t\t\t\t"tight'));
encode = strrep(encode, '"inward', sprintf('\t\t\t\t"inward'));
encode = strrep(encode, '"coordsys', sprintf('\t\t\t"coordsys'));
encode = strrep(encode, '"atlas', sprintf('\t\t"atlas'));
encode = strrep(encode, '"filepath', sprintf('\t\t\t"filepath'));

encode = strrep(encode, '"inputcoord', sprintf('\t\t\t"inputcoord'));
encode = strrep(encode, '"checkMRI', sprintf('\t\t"checkMRI'));
encode = strrep(encode, '"slices', sprintf('\t\t\t"slices'));
encode = strrep(encode, '"nslices', sprintf('\t\t\t"nslices'));

encode = strrep(encode, '"subj', sprintf('\t\t"subj'));
encode = strrep(encode, '"warpmni', sprintf('\t\t\t\t"warpmni'));
encode = strrep(encode, '"nonlinear', sprintf('\t\t\t\t"nonlinear'));
encode = strrep(encode, '"unit"', sprintf('\t\t\t\t"unit"'));

encode = strrep(encode, '"leadfield', sprintf('\t\t"leadfield'));
encode = strrep(encode, '"normalize', sprintf('\t\t\t"normalize'));
encode = strrep(encode, '"timeDomain', sprintf('\t\t"timeDomain'));

encode = strrep(encode, '"covariance', sprintf('\t\t\t"covariance'));
encode = strrep(encode, '"vartrl', sprintf('\t\t\t"vartrl'));
encode = strrep(encode, '"projectmom', sprintf('\t\t\t"projectmom'));
encode = strrep(encode, '"options', sprintf('\t\t"options'));
encode = strrep(encode, '"keep', sprintf('\t\t\t"keep'));
encode = strrep(encode, '"rawtrial', sprintf('\t\t\t"rawtrial'));
encode = strrep(encode, '"rep_timeseries', sprintf('\t\t"rep_timeseries'));

%% section 7: connectivity
encode = strrep(encode, '"connectivity', sprintf('\t"connectivity'));
encode = strrep(encode, '"filt_freqs', sprintf('\t\t"filt_freqs'));
encode = strrep(encode, '"collapse_band', sprintf('\t\t"collapse_band'));

%% dealing with special cases
% method: first two methods are indentend 3 times, last two are twice
% locations = strfind(encode, '"method"');
% newStr = regexprep('"', '"', '\t"');

%% save encode to correct location under correct name
% Write to a json file
fid = fopen(paths.mainconf, 'w');
fprintf(fid, '%s', encode);
fclose(fid);


















