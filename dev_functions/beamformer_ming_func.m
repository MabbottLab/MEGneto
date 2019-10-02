function [source,source_conn,parc_conn] = beamformer_ming_func(paths)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MEG PROCESSING SCRIPT - Beamforming and connectivity only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add paths to toolboxes and other dependencies

% addpath for Fieldtip toolbox

% addpath('/data2/dti_meg_ST/MEG/fieldtrip-20181231');
% ft_defaults
% 
% % addpath for analysis pipeline
% addpath(genpath('/data2/dti_meg_ST/MEG/MEGneto_2/'));
% 
% paths.base = '/data2/dti_meg_ST/MEG/rs_controls';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%
% fullPath = which('ft_preprocessing.m');
% [pathstr,~,~] = fileparts(fullPath);
% ftpath   = pathstr; % this is the path to FieldTrip

config = load_config(paths, paths.name);
config = config.config;
step = 'fcp2';
subj_ds = load_participants(paths,step);
pids = readtable(paths.all_subj_pids);
[subj_match, failure] = ds_pid_match(paths,step);
ssSubjPath = @(x) paths.(subj_match.pid{x});
rawdata = cellfun(@(x) [paths.rawdata '/' x], subj_match.ds, 'UniformOutput', false);
% info.ID = cellfun(@(x) [paths.anout, '/' x],subj_match.pid,"UniformOutput",false);
info.mri = cellfun( @(x) glob([paths.rawmri '/' x '*']),subj_match.pid);
info.ID = subj_match.pid;

mriformat = '*_V2.mri'; % can use any MRI format that includes fiducial coordinates (try CTFv2: '*_V2.mri')
ftpath = glob([find_megne2() '/external/fieldtrip*']);
ftpath = ftpath{1}(1:end-1);
% % Path to subject MRI file
% [LIST, ISDIR] = glob([paths.rawmri,'/' mriformat]);
% LIST(~ISDIR);
% %p.subj.subj_mri =  LIST;
% info.mri = LIST;
% 
% LIST = glob([paths.anout,'/*/']);
% temp = LIST;
% for i = 1:length(temp)
%     temp2 = strsplit(temp{i},'/');
%     ID(i) = temp2(end-1);
% end
% info.ID = ID;
% idx = ismember(info.ID, {'groupANALYSIS'});
% info.ID = info.ID(~idx)';
% 
% %numSub = length(data.ID)
% numSub = 1; % for testing 
%%
numSub = length(info.ID);
for i = 1:numSub
	inddata.mri = info.mri{i};
	inddata.ID = info.ID{i};
    fprintf('-------Running Participant : %s\n', inddata.ID);
	inddata.segmentationfigure = [ssSubjPath(i),'/',inddata.ID,'_segmentfigure.png'];
	inddata.sourcefigure = [ssSubjPath(i),'/',inddata.ID,'_sourcefigure.png'];
	%data = load preprocessind ica

	load([ssSubjPath(i),'/ft_meg_data_cfg.mat'])

	inddata.process = 'csd'; %Tlock (time-locked data 'old') or csd (cross spectral density) 
	% 'coh' = coherenace ; 'wpli_debiased' = debiased weighted phase lag index; 'granger' = granger causality; 'amplcorr' = amplitude correlation
	inddata.connectivitymethod = config.beamforming.connectivity.method;
    inddata.freq = [8 12]; %aplha [8 12] beta [13 29] theta [4 7]
    [source,source_conn,parc_conn] = fcp_3_beamforming_sourcegrid(inddata, ftpath, data, paths);
    disp('CONGRATUFUCKINGLATIONS YOU FINISHED ONE!!!!')
end

%%
