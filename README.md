# Free-Viewing MEGneto

This functional connectivity pipeline (fcp) is an add-on to the existing MEGneto fcp and is used for analysing Free-Viewing Data from the OIRM study. This Free-Viewing branch is intended to be used in tandem with the original MEGneto pipeline (i.e., both fcp's are required to conduct the analysis).

This pipeline was first developed and used throughout the months of June-August 2021. It was created for analyzing changes in neural oscillatory power in response to visual stimuli in the OIRM study. 

This pipeline is built on MATLAB using the FieldTrip toolbox. Developed @ SickKids Research Institute, Toronto, Canada. 

Please note that for extensive detail on how each step of the pipeline works, the user should navigate over to the main MEGneto repository's documentation, as the general idea of most steps (fcp_1-fcp_4 and fcp_5) remains consistent with the original pipeline.

- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [How to Use](#how-to-use)
   1. Extracting Clip Times
   2. Initial Setup
   3. JSON Config Setup
   4. Epoching
   5. Preprocessing
   6. ICA Checkpoint
   7. Channel Repair
   8. Beamforming
   9. Re-inserting Trials
   10. Frequency Analysis
   11. Statistical Analysis
- [Credits](#credits)

## System Requirements

* MATLAB
* FieldTrip Toolbox
* Machine with enough RAM

Note that, depending on available RAM on your system, the pipeline may crash during [beamforming](#beamforming) if your MEG data is not adequately downsampled or if you have requested too many virtual sources to be reconstructed (e.g., dipole grid resolution is too high). 

## Installation Guide

Download the repo through the Github website or use git in the command line to clone it on your machine. 

## How to Use

A template "main" function is provided under `templates/main_template.m` which guides the user through the pipeline steps. You should begin by making a copy of this file and renaming it. A unique main file should be created for each of your analyses, as it can serve as a record of what settings you used. 

A few important notes to remember before running the pipeline are:
1. The functions associated with the steps laid out below are found in the top-level MEGneto folder. Any related functions listed below are found in subfolders of the repo (e.g., the `functions` folder). 
2. The naming convention of your MRI files (which must have a .mri extension) is as follows. These file names should not have more than one underscore or period (i.e., the only period should be the file's extension `.mri`). If there is an underscore, the typical naming convention is `PID_version.mri`. 
3. Participant IDs follow the structure of study name prepended to participant number. For example, OIRM01 would mean that OIRM is the study name and 01 is the first participant. Other examples (for different studies) include ST05 and MEG04. The participant IDs do not require changing - they should be left as they already are. 
4. The pipeline can only process one task and condition at a time. If multiple tasks/conditions are fed in, there will be one set of .ds files for participants for task/condition 1 and one for task/condition 2, meaning there will be multiple .ds files for one participant. The pipeline is not equipped to handle this. If you have multiple tasks/conditions you wish to analyze, please do one at a time.

### Extracting Clip Times

Please navigate to the Extrapolating_clipTimes folder which contains two files (generate_markerTimes.m and read_data.m). 

The read_data.m file is the script for a function which is used in the generate_markerTimes.m script. 

The generate_markerTimes.m script is used to generate a .mat file containing the clip marks (also known as clip samples or clip times), the movie order for each run of MEG data, and the clippet marks (or samples/times). A pseudocode run-down of how this script operates can be found at the top of the script when it is opened. To run this file, you will need the OIRM dataset containing MEG data and .psy files (these files indicate the movie order for a given run of MEG data). You will also need a file that indicates the frame number where a clip/clippet change occured for each movie. This file should be provided to you and is titled main_cliptimes.mat or clipTimes.mat.

Note: keep in mind for conversion purposes that each frame is equal to 20 samples. 

Output: a *.mat file containing the clip/clippet times (the file's title should be something similiar to clipMarkers_allPpts.mat).

### Initial Setup

After making a copy of the main template and renaming it, open it and:
* Fill in the relevant folder paths and analysis name (lines 16-25)
* Add the Free-Viewing MEGneto pipeline, MEGneto and FieldTrip folders to the path, so MATLAB can find those functions (lines 29-31)
* Add a few additional folders to the path that are required to run the pipeline. This includes the path to the folder where the output from functions used to extrapolate clip times (line 35), the path to the folder where the feature vector which discriminates categories of visual content is located (line 36), and the path to folder where all OIRM MRI data is contained (line 37)
* Assuming this is your first run, skip the paths variable reloading at line 32
* Proceed to the section labelled "%%  fcp_0: setup (megne2setup)"

Please note that, due to the nature of the OIRM data set, it is advisable to conduct two separate analyses of the data - one for the "Run 1" data files and one for the "Run 2" data files. These two runs are combined after the beamforming step (freeviewing_fcp_4_beamforming.m).

`FREEVIEWING_MEGNE2SETUP.m` will create the folder structure for your analysis (e.g., config and analysis folders), and create empty setup files in the config directory (e.g., an empty config file to be filled in with analysis parameters, empty CSV files for each pipeline step that you will use to indicate which participants to analyze). If you already have a config file with your preferred parameters, replace the empty config file with that one. 

Note that this step will fail if your MEG and MRI data are not setup properly, namely:
* The `rawdata_path` string should be the path to the folder that contains all *.ds folders for all participants
* The `mri_path` string should be the path to the folder that contains all MRI files with the naming convention [PID].mri

Output: Struct called `paths` with all filepath definitions.

See also: 
- `freeviewing_path_generation.m` to generate path locations
- `path_check.m` to check that all paths are properly initialized

### JSON Config Setup

Prior to running the first step of the pipeline, the user must ensure that the JSON config file is populated with their desired parameters. `interactive_JSON_config.m` will prompt users to fill this JSON config file through an interactive graphical user interface (GUI). The user is repsonsible for filling in each field and sample inputs are presented to the user to demonstrate each field's format (note: the user can leave the sample input as is, if they wish to use that value for their analysis).

For the OIRM data set analysis conducted in the summer of 2021, the following fields were altered:
* email
* sampling rate (changd to 600Hz)
* the function field under taskFunc (changed to @freeviewingTaskTrialFun)
* the include field under tasktrialdef (changed to clipChange)
* the correct field under tasktrialdef (changed to clipChange)
* the template grid resolution (changed to 0.65)
* the atlas file path (changed to mmp)
* the first first in frequency analysis (changed to 1)
* the time window of interest in frequency analysis (changed to -1:0.033:3)
* the baseline in frequency analysis (changed to [-1, -0.5])
* the ROIs in frequency analysis 

For more detail on the meaning of each parameter in the JSON Config please see the main MEGneto repository's documentation.

### Epoching

`FREEVIEWING_FCP_1_TASKEPOCHING.m` will epoch MEG data into trials depending on the desired marker, detect trials with excessive head motion, muscle/jump artifacts, and bad channels. However, the epoching only rejects trials for excessive head motion and muscle/jump artifacts. Bad channels are detected and recorded, but repaired later on in the pipeline, after the ICA process at the final stage of preprocessing.

Within the task epoching step, the freeviewingTaskTrialFun is called. This function loads the output from the clip marker times extraction process (clipMarkerTimes.mat) and uses this file to epoch the data into trials which are marked by clip and clippet changes (lines 72-109). The function also ensures that all movie clips follow a consistent order for all particiapnts - that is, all run 1 data is ordered as movies 1-5 and all run 2 data is ordered as movies 6-10 (line 134). 

The task epoching step also identifies and outputs which trials are rejected. This may be important if the user wishes to average across trials later on in the pipeline and is thus required to re-insert blank columns at the appropriate indices for the missing trials. 

Note: if there are participants who do not have a matching MRI file, this step will not run until you: a) find the missing MRI and put it in the MRI folder, or b) remove their entry from the `subj_fcp1.csv`. 

Output: A struct with output file names, and for each subject: the number of trials per subject, trials marked with head motion, trials marked with noise, number of removed trials, indices of the removed trials, names of bad channels.

Notes:
- Ensure that subj_fcp1.csv is populated with the subject IDs of included participants.
- Prior to running this step, all desired parameters should be defined in the JSON config file. The user should double-check that the JSON config file is populated appropriately, especially if a template JSON was copied over. Information on the meaning of each parameter in the JSON config file can be found in the [Config Params Guide](https://github.com/dunjamatic/MEGneto/blob/configParams/ConfigParams.md).
- To see which lines differ from the fcp_1_taskepoching.m file in the MEGneto repository, please open the freeviewing_fcp_1_taskepoching.m file and find this information at the top of the script. 

### Preprocessing

`FREEVIEWING_FCP_2_PREPROCESSINGICA.m` will prepare epoched data for ICA preprocessing by downsampling and filtering with 3rd order gradients (derived from measurements taken by gradiometers). If indicated in the config JSON file, ICA will be carried out and the ICA components will be saved. The pipeline will downsample to whatever frequency the user specified in the config JSON.

Output: A struct with file names for the configuration of the preprocessed data, the data noise correlation matrix, and the ICA components.

Notes:
- Prior to running the function, ensure that `subj_fcp2.csv` is populated with the subject IDs of participants you want to include after checking over initial results.
- Outputs from fcp_1 will be loaded in at the start of this step. Additionally, a logging file will be set up to keep track of progress and the pipeline will check for matching MEG/MRI data. 
- Check participants who had excessive head motion or excessive numbers of bad channels.
- Need to remove bad channels from ica - if not you will get complex numbers. Because during repair channels procedure bad channels are repaired according to neighbours, thus the new ones are not unique (no independent components).
- To see which lines differ from the fcp_2_preprocessingica.m file in the MEGneto repository, please open the freeviewing_fcp_2_preprocessingica.m file and find this information at the top of the script. 

### ICA Checkpoint

`FREEVIEWING_FCP_2_5_CHECKPOINT.m` is an interactive session that guides the user through inspection of ICA components to identify components associated with artifacts such as heartbeats, blinks, etc. After inspection, the pipeline backprojects ICA components to remove the signal corresponding with the bad ICA components. For help on identifying components containing artifacts, see [ICA Inspection Guide](https://github.com/MabbottLab/MEGneto/blob/master/docs/ICA%20Inspection%20Guide%20v1.0.pdf).

Output: A struct with file names for the configuration of the preprocessed data, the data noise correlation matrix, and the ICA components and the bad components specified by the user. 

Notes:
- Prior to running the function, ensure that subj_fcp2_5.csv is populated with the subject IDs of participants you want to include.
- To see which lines differ from the fcp_2_5_checkpoint.m file in the MEGneto repository, please open the freeviewing_fcp_2_5_checkpoint.m file and find this information at the top of the script. 

### Channel Repair

`FREEVIEWING_FCP_3_CHANNELREPAIR.m` repairs bad channels detected from freeiviewing_fcp_1, but we held off on removing until the data had been ICA-cleaned. The channels are repaired by replacing them with some combination of neighbouring channels (default is 'weighted' average, other options include 'average', 'spline', or 'slap').

Output: *.mat file of fully cleaned data (i.e., removed head motion/muscle and jump artifacts, 3rd order gradients, ICA cleaned data, and repaired bad channels). 

Notes:
- Prior to running the function, ensure that subj_fcp3.csv is populated with the subject IDs of participants you want to include.
- Outputs from fcp2 are loaded at the start of this step and a logging file will be set up to keep track of progress. Further, the pipeline will check for matching MEG/MRI data.
- The output data of this step is in the sensor space MEG data (fully processed).
- To see which lines differ from the fcp_3_channelrepair.m file in the MEGneto repository, please open the freeviewing_fcp_3_channelrepair.m file and find this information at the top of the script. 

### Beamforming

`FCP_4_BEAMFORMING.m` maps functional data onto the source model and interpolates to an atlas. Here, be careful about conversions between mm and cm units in MEG and MRI data. A T1 template head model is loaded in (to normalize all participant head models), and is segmented to set the boundaries of which dipoles are actually located in the brain (ie. removes the skull). Once the head model is properly prepared to create the volume conduction model (which specifies how currents are generated by sources in the brain) a source model is prepared. The source model is a 3D grid/cortical sheet that specifies the set of positions of electric dipole currents that are considered in source reconstruction. 

For each participant, their anatomical MRI data is loaded in and fiducials are identified. The MRI data is segmented and prepared as a head model. The participant’s fully processed MEG data is loaded and resampled to the sampling rate specified in the config JSON file. A subject specific source model is prepared using the T1 template model as a control across participants. A leadfield is computed to be used for efficient inverse modelling. 

Finally, beamforming is performed. A covariance matrix is computed to inform us how related sensor signals are to each other, and sensor weights are computed using the gradiometers (sensor position), leadfield and head model. This information is used for beamforming. The participant’s trials are projected through a spatial filter (which transforms data from sensor to source space) and the resulting virtual sources are projected to their dominant orientation. The output at this point is a matrix with the timeseries for each source, trial, and participant. Next, an atlas is loaded into the pipeline and the source data is interpolated onto that atlas. In that interpolation process, all reconstructed sources that fall into regions as parcellated by the atlas are averaged or principal component analysis ([PCA](https://www.mathworks.com/help/stats/pca.html)) is performed and the first principal component is extracted, generating a representative time series for the region. 

Output: For each run of participant data, two matrices (atlas_beamforming_results_part1.mat and atlas_beamforming_results_part2.mat) than contain reconstructed timeseries for each trial and region of interest across both runs of MEG data. Thus, you can isolate the timeseries for a certain participant and trial etc.

Notes:
- Prior to running the function, ensure that subj_fcp4.csv is populated with the subject IDs of participants you want to include.
- At the start of the function, outputs from fcp_3 will be loaded in and a logging file will be set up to keep track of progress. Also, the pipeline will check for matching MEG/MRI data.
- To see which lines differ from the fcp_4_beamforming.m file in the MEGneto repository, please open the freeviewing_fcp_4_beamforming.m file and find this information at the top of the script. 

### Re-inserting Trials
`FREEVIEWING_FCP_4_5_REINSERTINGTRIALS.m` stacks the output of the beamforming step for each participant (i.e., it combines all beamforming from both runs of a participant's MEG data) and inserts blank columns for trials that were rejected in the task epoching step. 

Note: Please open the script of this code to see more details on this step.

Output: A matrix (with each run of a participant's data stacked together) with reconstructed timeseries for each trial and region of interest across both runs of MEG data for each participant. The output is stored in participant-specific folders under the name "revamped_beamforming_results.mat". 

### Frequency Analysis
`FREEVIEWING_FCP_5_FREQANALYSIS.m` uses spectral analysis on time-frequency representations of data to test hypotheses based on spectral power. The virtual sensor data from the reinserting trials step is loaded in and frequency analysis is performed on sliding timewindows of the data. Thus, for each subject and each interpolated atlas region, a power spectrum is calculated and corrected to a baseline to control for general/random spikes in power. 
For the free-viewing data analysis conducted in summer 2021, this function removes the previously inserted NaN trials (lines 117-128). The user is required to load in a feature vector that indicates which trials contain which category of visual content. Currently, this function is set up to load in a faces versus scenes feature vector (lines 83-85), and the user must specify in lines 135-137 which category of visual content is being analyzed (i.e., 0 for scenes, and 1 for faces) as the function outputs one result per category of visual content. If the user does not wish to separate by category of visual content they should specify that all trials should be included in the analysis (line 137).

Note: Frequency analysis requires the user to load in a feature vector (lines 

Output: A 4-D matrix containing power spectrum data. Matrix dimensions are [participants] x [regions] x [frequency] x [time]. Users can plot a power spectrum (frequency by time) for a specific region of a given participant’s data by plotting a slice of the matrix's first (participant) dimension (e.g., `imagesc(squeeze(powspctrm(participant_number, :, :, :)))`. 

Notes:
- Prior to running the function, ensure that subj_fcp5.csv is populated with the subject IDs of participants you want to include.
- At the start of this step, outputs from fcp_4 will be loaded in and a logging file will be set up to keep track of progress. Also, the pipeline will check for matching MEG/MRI data.
- To see which lines differ from the fcp_5_freqanalysis.m file in the MEGneto repository, please open the freeviewing_fcp_5_freqanalysis.m file and find this information at the top of the script. 

### Statistical Analysis
`FREEVIEWING_FCP_5_5_ANALYZEFREQANALYSIS.m` is used to perform statistical analysis and create visual representation of results for each aim and hypothesis of the OIRM Faces v Scenes Summer 2021 research project. This script is meant to be run in sections, not as a function, and the user should refer to the script for detailed instructions on each section of the script.

Note: Please open the script of this code to see more details on this step.

## Credits

**This list is incomplete!**

- Pre-2016: Marc Lalancette
- March 2016: Simeon Wong, Anne Keller
- November 2016: Sonya Bells
- June 2019: Ming Scott
- October 2019: Julie Tseng 
- August 2021: Dunja Matic
