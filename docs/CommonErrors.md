# A guide on commonly encountered errors and how to tackle them

- [Improper JSON config file set up](#improper-setup)
- [Subject .CSV file population](#populate-csv)
- [Subjects missing matching MRI and MEG data](#mismatch)
- [Incorrect file extension for MRI data (mri data is not .mri extension](#mri-filetype)
- [Incorrect interactive JSON config file population](#interactive-config)
- [Issues with file permissions on x11 forwarding](#x11-file-permissions)

## Improper JSON config file set up

### **Error occurrence #1:** During fcp_1, a warning pops up that says “not enough markers found”.

![](https://github.com/MabbottLab/MEGneto/blob/master/images/not_enough_markers.PNG)

Why did this occur?

In the JSON config file set up the user must specify marker names that are specific to
their data/analysis. Without this, when the data is read, the pipeline has no information
as to which events to epoch around (they are not found unless their names are
specified by the user).

How to debug:
1. Navigate to the JSON config file.
2. Find the `task.trialdef.markers` field and ensure these markers match the markers specific to your analysis. You should have a marker file that is specific to your data with which you can compare these fields.

Note: For more information on the JSON config file please navigate to the `ConfigParams.md` file for a description of each field.

### Error occurrence #2: Error in load config, parsing config, or recursive json struct.

![](https://github.com/MabbottLab/MEGneto/blob/master/images/JSON_population.PNG)

Why did this occur?

The user must populate a JSON config file which specifies various parameters, specific
to the user’s analysis/data, that the pipeline relies on for extracting information.
How to debug:

1. Navigate to the JSON config file and fill in the fields manually OR in the main
template file, run the `interactive_JSON_config` function which guides the
user through populating the config file.

Note: For more information on the JSON config file please navigate to the
ConfigParams.md file for a description of each field.
Forgetting to populate the subject .CSV files


## Subject .CSV file population

### Error occurrence #1: Error using load participants

![](https://github.com/MabbottLab/MEGneto/blob/master/images/CSV_population.PNG)

Why did this occur?

The user must populate the specific step’s comma separated value (csv) file with participants
such that the pipeline can perform the analysis on those participants. If that csv file is empty, the
pipeline has nothing to perform an analysis on/no data to work with.

How to debug:
1. Navigate to the step’s csv file (`subj_fcpX.csv`, where X is the number of the step you
are on in the pipeline) and populate the file with participants OR uncomment the lines
above the call to fcp1 to auto populate the csv file. If you wish to auto populate
subsequent steps, copy the code above the call to fcp1 and change the `subj_fcp1` is
`fid = fopen(paths.subj_fcp1, 'w')` to match the step number that you are on.

    a. Note that auto population is fine to do for every step, however the user must
ensure they manually delete any participants they wish to exclude (due to an
excess of bad channels, bad trials, etc.).

### Error occurrence #2: Using auto-population beyond fcp1

Why did this occur?

Prior to fcp1, it is common for the user to user the pipeline's code to auto-populate the subj_fcp1.csv
file. However, users should not repeat this for subsequent steps since in most analyses, some participant's 
will be removed after several steps in the pipeline. 

Instead of using the auto-population code,
user's should navigate to the previous fcpX step's subject CSV file and copy the participant list to the 
current fcpX steps subject CSV file, then remove any more participants that they need to.

How to debug:
1. Do not use the auto-population code for steps other than fcp1. Instead, do the steps listed below.
2. Navigate to the previous step’s csv file (`subj_fcpX.csv`, where X is the number of the step you
just finished in the pipeline) and copy the list of participants.
3. Navigate to the current step’s csv file (`subj_fcpX.csv`, where X is the number of the step you
are about to run in the pipeline) and paste the list of participants. Remove any participant's you 
do not wish to include in the analysis.

## Subjects missing matching MRI and MEG data

### Error occurrence : Error using ds_pid_match which finds matching MEG and MRI data.

![](https://github.com/MabbottLab/MEGneto/blob/master/images/MEG-MRI_mismatch.PNG)

Why did this occur?

MEGneto requires MRI data to have the .mri extension. If a participant’s MRI data that does not
have this extension, the pipeline will error out.

How to debug:
1. Ensure all MRI files end in .mri. If they do not, check if there is existing .mri data for
these participants elsewhere and add it to the MRI folder. Else, delete the participants
that do not have .mri data.

## Incorrect MRI file naming convention

### Error occurrence: Error occurs in megne2setup’s path generation.

![](https://github.com/MabbottLab/MEGneto/blob/master/images/naming_convention.PNG)

Why did this occur?

Participant MRI data follows a naming convention. If the naming convention is violated, the
pipeline errors out as it gets confused.

How to debug:
1. Navigate to the folder which contains the MRI data and alter the file names to match the
following convention
a. MRI files should only contain one underscore after the last slash that designates
the file path (e.g., `/xxx/xxx/xxx/ST01_V1.mri` is acceptable but
`/xxx/xxx/xxx/ST01_V1_V2.mri` is not). Note that underscores are not
required (e.g., `/xxx/xxx/xxx/ST01.mri` is acceptable).

## Incorrect interactive JSON config file population

### Error occurrence: Incorrect format inputted to interactive JSON config 

Why did this occur?
The interactive config file presents a template/sample input for each field. The user
is welcome to change any of these values, however they must follow the format
presented in the interactive JSON config, else errors will occur.

How to debug:
1. Ensure you are following the exact format of inputs presented in the interactive config.Common mistakes include putting spaces between characters that shouldn't have spaces (e.g., freqanalaysis.ROIs is inputted as [1,2,3];[4];[5,6] with no spaces between any of the characters).


## Issues with file permissions on x11 forwarding

### Error occurrence: Cannot save output of interactive JSON config
![](https://github.com/MabbottLab/MEGneto/blob/master/images/x11_error.png)

Why did this occur?
The JSON file does not have write permissions, so it cannot save the changes the user is trying to make to it.

How to debug:
1. Using powershell, navigate to the folder which contains the JSON config file (it is located within the config folder). As a verification that the issue is indeed a lack of writing permission, type `ls -l` and you will see the file permissions for each file in this folder, as pictured below.
![](https://github.com/MabbottLab/MEGneto/blob/master/images/x11_permission_verification.png)

The information displayed in the leftmost column on the image represents the permissions in clusters, as described in the image below.
![](https://github.com/MabbottLab/MEGneto/blob/master/images/permission_description.png)

As seen in the first image, the JSON file (`OIRMFake.json`, highlighted in green), does not have a w in the owner cluster, indicating that there are no writing permissions to this file. To remedy this, simply type the following command, while remaining in this folder: `chmod +w OIRMFake.json` (where OIRMFake.json would be the name of your JSON file).
