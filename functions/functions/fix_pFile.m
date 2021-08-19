clear

%%
load('/home/sonya/Data_exercise/T3_I/groupANALYSIS/pINFO.mat','p'); % load p-structure

%%
p.paths.output_dir = @(id) [p.subj.subj_dir{id},'ANALYSIS/'];

%%% INPUTS / OUTPUTS %%% (update anonymous function's workspace if any edits to subject paths made)
p.paths.epoching_info = @(id,ext) [p.paths.output_dir(id),'epoching_',num2str(p.thr),'mmHM.',ext];
p.paths.trial_cfg = @(id) [p.paths.output_dir(id), 'ft_meg_trl_cfg.mat']; 
p.paths.grad_cfg =   @(id) [p.paths.output_dir(id), 'ft_meg_grad_cfg.mat'];
p.paths.subj_epochInfo  =  [p.paths.group_output_dir,'subj_epoching_info.mat']; 
p.paths.group_rmBadChan =         [p.paths.group_output_dir, 'group_rmBadChan.mat'];

p.paths.preprocessedData_cfg  =  @(id) [p.paths.output_dir(id), 'ft_meg_data_cfg.mat']; 
p.paths.preprocessedData_grad  =  @(id) [p.paths.output_dir(id), 'ft_meg_data_grad.mat']; 
p.paths.ICAcomp_cfg           =  [p.paths.group_output_dir,'ft_icacomp.mat']; % save ica components
p.paths.subj_icafigure        =  @(id,ext) [p.paths.output_dir(id), p.subj.ID{id}, '_icafigure.', ext];

% Final output for virtual sensors (116 AAL region timeseries)
p.paths.subj_116vsoutput =          @(id,postfix) [p.paths.output_dir(id), p.subj.ID{id},'_116VS_', postfix, '.mat'];
% Final output for virtual sensors (all source timeseries)
p.paths.subj_vsoutput =             @(id,postfix) [p.paths.output_dir(id), p.subj.ID{id}, '_VS_', postfix, '.mat'];
% Final output for virtual sensor coordinates (all source timeseries)
p.paths.subj_source_subjcoords =    @(id,postfix) [p.paths.output_dir(id), p.subj.ID{id}, '_VS_', postfix, '_subjcoords.mat'];
% Final output for virtual sensor labels (all source timeseries)
p.paths.subj_source_labels =        @(id,postfix,ext) [p.paths.output_dir(id)', p.subj.ID{id}, '_VS_', postfix, '_labels.', ext];
% Output for head model (matlab structure)
p.paths.subj_headmodel =            @(id,postfix) [p.paths.output_dir(id), p.subj.ID{id}, '_headmodel_', postfix, '.mat'];
% Output for source model (matlab structure)
p.paths.subj_sourcemodel =          @(id,postfix) [p.paths.output_dir(id), p.subj.ID{id}, '_sourcemodel_', postfix, '.mat'];
% Output for full source analysis (matlab structure)
p.paths.subj_sourceanalysis =       @(id,postfix) [p.paths.output_dir(id), p.subj.ID{id}, '_sourceanalysis_', postfix, '.mat'];
% Source Figure
p.paths.subj_sourcefigure =         @(id,postfix,ext) [p.paths.output_dir(id), p.subj.ID{id}, '_sourcefigure_', postfix, '.', ext];
% Segmentation Figure
p.paths.subj_segmentationfigure =   @(id,postfix) [p.paths.output_dir(id), p.subj.ID{id}, '_segmentationfigure_', postfix, '.png'];

%%% FIDUCIAL PARSE FUNCTION
% strip non-numeric / dash / period / space characters
% convert to numeric array
% reshape to column vector
p.fiducial.parse = @(fidstr) reshape(str2num(fidstr(fidstr >= 48 & fidstr <= 57 | fidstr == 45 | fidstr == 46 | fidstr == 32)), [], 1); %#ok<ST2NM>

p.paths.connmatNetwork  = @(id,postfix,network) [p.paths.output_dir(id), p.subj.ID{id}, '_', postfix, '_', network,'.mat'];
        
        
p.paths.subj_116vsoutput =  @(id,postfix) [p.paths.output_dir(id), p.subj.ID{id},'_116VS_', postfix, '.mat'];

%%
% save p-structure
disp('Saving...');
save(p.paths.p_strct,'p','-mat','-v7');
disp('Done.');
