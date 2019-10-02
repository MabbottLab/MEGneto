%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MEG PROCESSING SCRIPT - Beamforming and connectivity only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add paths to toolboxes and other dependencies

% addpath for Fieldtip toolbox

addpath('/data2/dti_meg_ST/MEG/fieldtrip-20181231');
ft_defaults

% addpath for analysis pipeline
addpath(genpath('/data2/dti_meg_ST/MEG/MEGneto_2/'));

paths.base = '/data2/dti_meg_ST/MEG/rs_controls';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fullPath = which('ft_preprocessing.m');
[pathstr,~,~] = fileparts(fullPath);
ftpath   = pathstr; % this is the path to FieldTrip

paths.mri = [paths.base,'/*/*MRI*'];
mriformat = '*_V2.mri'; % can use any MRI format that includes fiducial coordinates (try CTFv2: '*_V2.mri')

% Path to subject MRI file
[LIST, ISDIR] = glob([paths.mri,'/' mriformat]);
LIST(~ISDIR);
%p.subj.subj_mri =  LIST;
info.mri = LIST;

[LIST, ISDIR] = glob([paths.base,'/*/']);
temp = LIST;
for i = 1:length(temp)
    temp2 = strsplit(temp{i},'/');
    ID(i) = temp2(end-1);
end
info.ID = ID;
idx = ismember(info.ID, {'groupANALYSIS'});
info.ID = info.ID(~idx)';

%numSub = length(data.ID)
numSub = 1; % for testing 
%%
 i = numSub;
for i = 1:numSub
	inddata.mri = info.mri{i};
	inddata.ID = info.ID{i};
    fprintf('-------Runing Participant : %s\n', inddata.ID);
	inddata.segmentationfigure = [paths.base,'/',inddata.ID,'/',inddata.ID,'_segmentfigure.png'];
	inddata.sourcefigure = [paths.base,'/',inddata.ID,'/ANALYSIS/',inddata.ID,'_sourcefigure.png'];
	%data = load preprocessind ica

	load([paths.base,'/',info.ID{i},'/ANALYSIS/ft_meg_data_cfg.mat'])

	inddata.process = 'csd'; %Tlock (time-locked data 'old') or csd (cross spectral density) 
	% 'coh' = coherenace ; 'wpli_debiased' = debiased weighted phase lag index; 'granger' = granger causality; 'amplcorr' = amplitude correlation
	inddata.connectivitymethod = 'coh';
    inddata.freq = [8 12]; %aplha [8 12] beta [13 29] theta [4 7]
    [source,source_conn,parc_conn] = fcp_3_beamforming_sourcegrid(inddata, ftpath, data);
end

%%
