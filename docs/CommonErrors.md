# A guide on commonly encountered errors and how to tackle them

- [Improper JSON config file set up](#improper-setup)
- [Forgetting to populate the subject .CSV files](#populate-csv)
- [Subjects missing matching MRI and MEG data](#mismatch)
- [Incorrect file extension for MRI data (mri data is not .mri extension](#mri-filetype)

## Improper JSON config file set up

### **Error occurrence #1:** During fcp_1, a warning pops up that says “not enough markers found”.

[@DUNJA: INSERT SCREENSHOT HERE]

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

[@DUNJA: INSERT SCREENSHOT HERE]

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

### Error occurrence : Error using load participants

[@DUNJA: INSERT SCREENSHOT HERE]

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

## Subjects missing matching MRI and MEG data

### Error occurrence : Error using ds_pid_match which finds matching MEG and MRI data.

Why did this occur?

MEGneto requires MRI data to have the .mri extension. If a participant’s MRI data that does not
have this extension, the pipeline will error out.

How to debug:
1. Ensure all MRI files end in .mri. If they do not, check if there is existing .mri data for
these participants elsewhere and add it to the MRI folder. Else, delete the participants
that do not have .mri data.

## Incorrect MRI file naming convention

### Error occurrence: Error occurs in megne2setup’s path generation.

[@DUNJA INSERT SCREENSHOT HERE]

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
