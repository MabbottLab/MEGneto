# MEGneto 3.0 

This functional connectivity pipeline (fcp) is built on MATLAB using the FieldTrip toolbox to analyze MEG data. Developed @ SickKids Research Institute, Toronto, Canada. See the PDF 'workflow' to get an overview of the pipeline.

- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [How to Use](#how-to-use)
   1. [Initial Setup](#initial-setup)
   2. [Epoching](#epoching)
   3. [Preprocessing](#preprocessing)
   4. [ICA Checkpoint](#ica-checkpoint)
   5. [Channel Repair](#channel-repair)
   6. [Beamforming](#beamforming)
   7. [Functional Connectivity](#functional-connectivity)
- [Credits](#credits)
- [On Downsampling](#on-downsampling)

## System Requirements

This pipeline is currently being developed with MATLAB R2019a in a Linux environment. Configuration is set using JSON files, inspected with a any basic text editor. The FieldTrip toolbox (also required) contains compatibility functions should you need older or newer versions of certain key functions. 

Note that, depending on available RAM on your system, the pipeline may crash during beamforming (fcp_3) if your MEG data is not adequately downsampled. (See the [On Downsampling](#on-downsampling) section for more on how to handle this.)

## Installation Guide

Download the repo or clone it on your machine in a sensible place. 

## How to Use

Before you run the pipeline, this repo must be fully visible in the path, as well as the top-level FieldTrip folder. You can do this by running the following lines:

``` MATLAB
addpath(genpath('/path/to/MEGneto'))
addpath('/path/to/FieldTrip') % note the lack of genpath here
ft_defaults; % allow fieldtrip to run setup
```

Note that the functions associated with the steps laid out below are found in the top-level MEGneto folder. Any related functions listed below are found in subfolders of the repo (e.g., the `functions` folder). Anything under development is, accordingly, under `dev_functions`. 

*DOCUMENTATION UNDER CONSTRUCTION - NEEDS A BIG UPDATE, I'M WORKING ON IT OKAY?? - JT, 2020-04-03*

*I'm reducing the amount of detail here to provide a higher level skeleton rather than play-by-play.*
*Especially since the code is now much more thoroughly documented. -JT, 2020-04-06*

### Initial Setup

`MEGNE2SETUP.m` will create a subfolder to `PROJECT_PATH` or the current working directory (if `PROJECT_PATH` is not provided) called ANALYSIS_NAME, create config and analysis subfolders within, and create unfilled `participants.txt` and `config.json` files in the config directory. If you already have a config file with your preferred parameters, replace the unfilled config file with that one. You should have a rawdata folder in your project directory with all the `*.ds` files, and a subfolder to rawdata called `MRIs` which contains all the MRIs in the form `[PID].mri`.

See also: `path_generation.m`, `path_check.m`

### Epoching

`FCP_1_TASKEPOCHING.m` will epoch MEG data into trials depending on the desired marker, detect trials with excessive head motion, muscle/jump artifacts, and bad channels. However, the epoching only rejects trials for excessive head motion and muscle/jump artifacts. Bad channel repair happens later, after the ICA process at the final stage of preprocessing.

Notes: 
- Ensure that `subj_fcp1.csv` is populated with the subject IDs of included participants. 
- Desired parameters should be defined in the JSON config file. User should double-check that the JSON config file is populated appropriately, especially if a template JSON was copied over. 

See also: `DS_PID_MATCH`, `WRITE_MATCH_IF_NOT_EMPTY`, `PLOTTRIGGERS`, `HEADMOTIONTOOL`, `DETECTBADCHANNELS`, `SEARCHTASKTRIALFUN`, `DETECTBADCHANNELS`

### 2.  Preprocessing

`FCP_2_PREPROCESSINGICA.m` will prepare epoched data for ICA preprocessing by downsampling and filtering with 3rd order gradients, then carry out the ICA omitting signal from bad channels. The pipeline will downsample to whatever frequency you've specified in the config JSON.  

Notes:
- Check participants who had excessive head motion or excessive numbers of bad channels. 
- Ensure that subj_fcp2.csv is populated with the subject IDs of participants you want to include after checking over initial results. 
- Need to remove bad channels from ica - if not you will get complex numbers. Because during repair channels procedure bad channels are repaired according to neighbours, thus the new ones are not unique (no independent components).

See also: `ft_channelselection.m`, `ft_selectdata.m`, `ft_resample.m`, `ft_componentanalysis.m`, 

### ICA Checkpoint

`FCP_2_5_CHECKPOINT.m` is an interactive session that guides user through inspection of ICA components to identify components associate with artifacts such as heartbeats, blinks. After inspection, backprojects ICA components to remove signal corresponding with bad ICA components. 

See also: `disp_ica_chans.m` (at bottom of script), `ft_rejectcomponent.m`, `ft_databrowser`

### Channel Repair

`FCP_3_CHECKPOINT.m` repairs bad channels detected from fcp_1, but we held off on removing until the data had been ICA-cleaned. 

See also: `ft_prepare_neighbours`, `ft_channelrepair`

### Beamforming

`FCP_4_BEAMFORMING.m` maps functional data onto the source model and interpolates to AAL atlas regions (currently excluding cerebellar regions). Here, be careful about conversions between mm and cm units in MEG and MRI data. 

See also: 
- `ft_read_mri` to import T1 template from spm8
- `ft_volumesegment.m` to segment anatomical MRI into T1 template specs
- `ft_prepare_headmodel.m` to construct a volume conduction model based on geometry of head, takes previous output as input
- `ft_convert_units.m` to convert volumes between mm and cm for CTF type
- `ft_read_atlas.m`
- `ft_volumelookup.m` to create binary mask; once applied, will isolate desired regions
- `ft_sourceplot.m` to visualize/check alignment between template and subject head models; save output image
- `ft_prepare_leadfield.m` to compute the lead field matrix
- `ft_sourceanalysis.m` for source reconstruction
- `ft_sourcedescriptives.m` to project to dominant orientation (largest eigenvector)
- `ft_sourceinterpolate.m` to interpolate functional data onto anatomical data using prev as input, subject MRI

### Functional Connectivity

`FCP_5_TASKCONNECTIVITY.m` estimates functional connectivity. Currently supports the following connectivity metrics: 
- 'plv' (phase locking value)
- 'pli' (phase lag index)
- 'wpli' (weighted phase lag index)
- 'wpli_debiased' (as above, but debiased)
- 'coh' (coherence)
The strings listed above should be specified in the config JSON file exactly as presented. 

See also: `ft_freqanalysis.m`, `ft_connectivityanalysis.m`, `ft_checkdata.m`

## Credits

**This list is incomplete!**

- March 2016: Simeon Wong, Anne Keller
- November 2016: Sonya Bells
- June 2019: Ming Scott
- October 2019: Julie Tseng 

## On Downsampling

You may run out of RAM during the beamforming step if your MEG data is not adequately downsampled. MEG data is typically gathered at 1200Hz or 600Hz. However, based on the [Nyquist-Shannon sampling theorem](https://en.wikipedia.org/wiki/Nyquist%E2%80%93Shannon_sampling_theorem), you can downsample to 2x the max frequency you want to recover. But what does that mean in practice?

Assume that the max frequency you want to analyze is 100Hz (e.g., the upper limit of the high gamma frequency band). Then, the Nyquist-Shannon theorem suggests that your data needs to be sampled at 100 x 2 = 200Hz to recover information properly from that 100Hz frequency. 

Based on this, you may choose to downsample your data from 1200Hz to 300Hz. For a 4-second epoch (e.g., -2s to +2s around a marker of interest), this means going from 4800 timepoints to 1200 timepoints for each of the *151* channels. That's 724,800 points to 181,200 - almost 600k less timepoints to crunch with no information loss. Way easier on the machine. 
