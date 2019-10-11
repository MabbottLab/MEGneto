# MEGneto: MEGne3 Branch
 
Pipeline on MATLAB with FieldTrip to process MEG data and run functional connectivity analyses. Developed @ SickKids, Toronto, Canada. 

See the PDF 'workflow' to get an overview of the pipeline.

## Modules

Important: this repo and FieldTrip top-level folder should be visible upfront.

### 0. Setup

   - *megne2setup.m*: creates project folders incl. output
      - Calls *ft_defaults.m*: FieldTrip initialization that adds proper subfolders
      - Checks input params (folder names, data existence)
      - Calls *path_generation.m*: setup folder/file names
      - Creates folders according to path_generation output
      - Calls *path_check.m*: verify successful folder/file creation

#### Critical fixes
- megne2setup.m, line 33: Add path to FieldTrip on system
   FieldTrip used to be held within this repo but removed for space efficiency.

#### Nice-to-haves
- megne2setup.m, line 36: Comb through ft_defaults.m to customize settings according to our needs (e.g., turn off track usage, specify using MATLAB toolboxes rather than compat
- megne2setup.m, line 63: change folder structure of output (adjust this checkpoint as to whether the results already exist)
- path_generation.m: specify output folder structure
- megne2setup.m: generally, clean up excessive data type conversion

### 1. Preprocessing | fcp_1_TaskEpoching.m

Detects and removes timewindows or full trials with excessive head motion, muscle/jump artifacts, detects noisy channels but does not remove. 

- Calls *load_config.m*: read in config settings
- Calls *load_participants.m*: grab *.ds file locations from subj_fcp1.csv (or throw error if empty)
- Calls *ds_pid_match.m*: match *.ds and *.mri subject IDs
- Check that there are matched pairs of MEG/MRI files
- Calls *write_match_if_not_empty.m*: record matched subject ID list, incl. if number of matches has changed OR on first pass, backs up previous list otherwise
- Initializes several output files (e.g., head motion PNG, subject epoching info in a *.mat file)
- Calls *save_to_json.m*: save output meta-data
- Calls *plotTriggers.m*: visualize events such as button press, fixation cross presentation; saves triggers
- Calls *ft_definetrial.m*: FieldTrip function that defines trials for further processing and analysis
   - Calls *searchTaskTrialFun.m*: custom trial function that discretizes continuous event lists into trials, identifies successful responses or other desired markers
- Calls *HeadMotionTool.m*: Identify bad trials with motion over specified threshold in config file; saves tool image output to visualize head motion
- Remove those trials and save output in JSON
- Specifies muscle artifact filtering parameters
- Calls *ft_artifact_muscle.m*: Returns windows corresponding to muscle artifacts
- Specifies jump artifact threshold
- Calls *ft_artifact_jump.m*: same idea as muscle ver
- Specifies whether complete or partial trial rejection should occur
- Calls *ft_rejectartifact.m*: removes trials or trial segments corresponding to muscle and jump artifacts
- Assembles a Nx2 array of timewindows corresponding to both classes of artifacts
- Specifies bad channel detection parameters
- Calls *detectBadChannels.m*: setup to bad channel detection, then
   - Calls *detectBadChannels_Algorithm.m*: actual calculations relevant to noisy channel detection; returns 1 for noisy, 0 otherwise
- Save output to JSON file w/in ./analysis/PID/ folders

#### Critical fixes
- fcp_1_TaskEpoching.m: Add capability to specify what class of trigger (e.g., left correct or right correct?)
- PID *.ds and *.mri matching is messy right now with unstable folder structure

#### Nice-to-haves
- fcp_1_TaskEpoching.m: Input multiple triggers at once to epoch for, loop over
- fcp_1_TaskEpoching.m: Option to epoch before/after artifact removal, noisy channel removal
- JSON template
- fcp_1_TaskEpoching.m, line 133: I think the config settings are already found within the JSON config, then overwritten by this section.
- JSON config setting of 1 for rmBadChannels meaningless - currently runs it regardless of specification, and without any options to modify bad channel detection parameters

### 2.  Preprocessing | fcp_2_PreprocessingICA.m

Allows for manual removal of problematic channels, preparation and execution of ICA, user guided removal of ICA components that reflect heartbeat/blink artifacts. 

- Pull in cfg JSONs that contain epoch information
- Import and filter raw *.ds data
- Load gradiometer config file
- Load 3rd order gradients for noise reduction in CTF
- Specify whether bad channel removal is manual or not?
   - If manual, call *inputBADchannels.m*: interactive user-input
- Call *ft_prepare_neighbours.m*: Specify sensor removal parameters, prepare
- Call *ft_channelrepair*: replaces bad or missing channels with plain average of neighbouring sensors
- If no ICA:
   - Resample and return
- If yes ICA:
   - Grab list of bad channels
   - Call *ft_channelselection.m*: generate list of channels to incorporate
   - Call *ft_selectdata.m*: isolate data based on prior output
   - Specify ICA parameters (e.g., runica, fastica, etc. - currently set to runica without specifying PCA num components)
   - Call *ft_componentanalysis.m*: run ICA
   - Save output
   - If interactive: 
      - Display each component for review and record user input
   - If not itneractive:
      - Get list of components to be removed
      - Call *ft_rejectcomponent.m*
      - Save output

#### Critical fixes
- Add bad channels output from fcp_1 to be removed here?

#### Nice-to-haves
- Save filtered ft_preprocessing output to prevent repetitions
- Separate bad channel removal or channel repair from ICA
- Specify PCA to run before ICA and num components to reduce run time

### 3. Beamforming

Source reconstruction analysis.

- Load config information from JSON files, check IDs
- Call *ft_read_mri* to import T1 template from spm8
- Specify coordinate system in FT
- Call *ft_volumesegment.m*: segments anatomical MRI into T1 template specs
- Call *ft_prepare_headmodel.m*: constructs a volume conduction model based on geometry of head, takes prev as input
- Call *ft_convert_units.m*: convert volume to cm for CTF type
- Create figure w/ template head model, dipole grid
- Call *ft_read_atlas.m*, *ft_convert_units.m*: load atlas, convert units
- Call *ft_volumelookup.m*: create binary mask; once applied, will isolate desired regions
- Load subject's MRI, preprocessed MEG data
- Carry out similar procedure to template head model generation above
- Call *ft_sourceplot.m*: visualize/check alignment between template and subject head models; save output image
- Check alignment between subject head and source model; save output image
- Call *ft_prepare_leadfield.m*: compute the lead field matrix
- Source reconstruction:
   - Either 'Tlock' or 'csd'
   - Set parameters
   - Call *ft_sourceanalysis.m* + other steps specific to the type of source recon
- Call *ft_sourcedescriptives.m*: project to dominant orientation (largest eigenvector)
- Call *ft_sourceinterpolate.m*: interpolate functional data onto anatomical data using prev as input, subject MRI

#### Critical fixes
- Separate connectivity analysis into its own code

#### Nice-to-haves
- Ability to parallel process participants through analysis

### 3. Functional connectivity analysis

- Call *ft_connectivityanalysis.m*

## Historical contributors

**This list is incomplete!**

- March 2016: Simeon Wong, Anne Keller
- November 2016: Sonya Bells
- June 2019: Ming Scott
- October 2019: Julie Tseng 