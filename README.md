# MEGneto: MEGne3 Branch
 
Pipeline on MATLAB with FieldTrip to process MEG data and run functional connectivity analyses. Developed by the Mabbott Lab @ SickKids, Toronto, Canada. 

Important: this repo and FieldTrip top-level folder should be visible upfront.

## Modules

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

### 1. Preprocessing

   - *fcp_1_TaskEpoching.m*: epoching
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



#### Critical fixes
- fcp_1_TaskEpoching.m: Add capability to specify what class of trigger (e.g., left correct or right correct?)
- PID *.ds and *.mri matching is messy right now with unstable folder structure

#### Nice-to-haves
- fcp_1_TaskEpoching.m: Input multiple triggers at once to epoch for, loop over
- fcp_1_TaskEpoching.m: Option to epoch before/after artifact removal, noisy channel removal
- JSON template
- fcp_1_TaskEpoching.m, line 133: I think the config settings are already found within the JSON config, then overwritten by this section.
- JSON config setting of 1 for rmBadChannels meaningless - currently runs it regardless of specification, and without any options to modify bad channel detection parameters

### 2. Beamforming

### 3. Functional connectivity analysis

## Historical contributors

- March 2016: Simeon Wong, Anne Keller
- November 2016: Sonya Bells (November 2016)
- June 2019: Ming Scott
- October 2019: Julie Tseng 