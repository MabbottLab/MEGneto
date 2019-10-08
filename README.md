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

#### Would-be-nice
- megne2setup.m, line 36: Comb through ft_defaults.m to customize settings according to our needs (e.g., turn off track usage, specify using MATLAB toolboxes rather than compat
- megne2setup.m, line 63: change folder structure of output (adjust this checkpoint as to whether the results already exist)
- path_generation.m: specify output folder structure
- megne2setup.m: generally, clean up excessive data type conversion

### 1. Preprocessing

### 2. Artifact removal

### 3. Beamforming

### 4. Functional connectivity analysis


