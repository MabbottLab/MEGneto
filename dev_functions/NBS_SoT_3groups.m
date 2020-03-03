%% To make data matrix and design matrix for NBS exercise 
% modified from NBS_sam by Sonya - May 9, 2018
% modified by Liz - Sept 2018, Jan 2020

% load path to parpticipant data matrix files for participants
spath = '/mnt/sda/juanita/MEGneto/analysis/Right/analysis';

%% create data matrix - PLI adj. matrix (subjects within conditions)

% load pli matrix files for each participant in CONTROL 
% grouPLI{} = LIST
[LIST, ISDIR] = glob([spath,'/ST*/fcp_5_adjmat_wpli_deb.mat']);
% [LIST, ISDIR] = glob([spath,'/control/ST*/fcp_5_adjmat_wpli_deb.mat']); % here we are using debiased wPLI
LIST(~ISDIR);
groupPLI{1} = LIST;
% group1PLI = LIST; % group1PLI = controls

%load pli matrix files for each participant in SURGERY ONLY
% [LIST, ISDIR] = glob([spath,'/surg/ST*/fcp_5_adjmat_wpli_deb.mat']);
[LIST, ISDIR] = glob([spath,'/ST*/fcp_5_adjmat_wpli_deb.mat']);
LIST(~ISDIR);
groupPLI{2} = LIST; % group2PLI = surgery only

%load pli matrix files for each participant in RAD/CHEMO 
[LIST, ISDIR] = glob([spath,'/rad/ST*/fcp_5_adjmat_wpli_deb.mat']);
LIST(~ISDIR);
groupPLI{3} = LIST; % group3PLI = surgery + rad and/or chemo


%% initialize subjects in conditions - load PLI adj. matrix
% 
% numGroups = 3; % number of groups 
% gg = 1; % group 1 
% numCond = 1; % number of conditions 
matlinks = [];
numGroups = 3;
for gg = 1:numGroups
       fprintf(['\t---- Group ',num2str(gg),' ----\n']);
            for ss = 1:length(groupPLI{gg})
                this_adjmat = load(groupPLI{gg}{ss});
                matlinks = cat(3, matlinks, this_adjmat.adjmat(:,:,1,:));
            end
end

%% Save matrix by frequency band 

network.file_out= '/mnt/sda/juanita/MEGneto/analysis/left/analysis';
freq ={'theta','alpha','beta','Lgamma','Hgamma'};
nbs_datamat = matlinks;
% save NBS data matrix
for ff= 1:5 % range freq save
    data_matrix = nbs_datamat(:,:,:,ff);
    datamat_output = [network.file_out,'/NBS',freq{ff},'_Hz.mat'];
    save(datamat_output,'data_matrix','-mat');
    fprintf(['NBS Data Matrix saved to: ',datamat_output,'\n']);
end

    %datamat_output = [network.file_out,'/NBS_',num2str(network.filt_freqs(ff,1)),'-',num2str(network.filt_freqs(ff,2)),'Hz.mat'];

%% make design matrix for NBS
% modified from Julie Sato's NBS script (written by Simeon Wong)

numSubj = 1:length(groupPLI{ss});

design_matrix = zeros(size(nbs_datamat,3),3);
design_matrix(1:length(numSubj),1) = 1;
design_matrix(11:end,1) = 0; % change index based on where second group begins
design_matrix(11:end,2) = 1; 

% design_matrix = zeros(size(nbs_datamat,3),2);
% design_matrix(1:length(numSubj),1) = 0;
% design_matrix(18:end,1) = 1;
% design_matrix(1:length(numSubj),2) = -1;
% design_matrix(18:end,2) = 1;


% save 
designmat_output = [network.file_out,'/design_matrix.mat'];
save(designmat_output,'design_matrix','-mat');
fprintf(['\nDesign Matrix saved to: ',designmat_output,'\n\n']);

% below: extra script from fcp_5_NBS.mat in /connectivity_study/analysis_tools/MEGneto

% UI.contrast.ui=[1,-1];                % Numeric array specifying contrast
%                                             % Example:      Comparing two groups
%                                             %                    -and [1,-1]  [-1,1] - [G1 G2] - 'Group'Differences
%                                             % Example:      Comparing two 'time points' or conditions
%                                             %                    - [1,-1] and [-1,1] - [C1 C2] - 'Condition'Differences
%                                             % Example:      Comparing two 'time points' or conditions in two groups
%                                             %                    - [1,-1] and [-1,1] - [G1C2-G1C1 G2C2-G2C1] - 'Conditions in Group'Differences
% network.difference = 'Condition';     % Difference Interested in 'Group' | 'Condition' | 'Conditions in Group'


% if numCond > 1
%     % paired groups
%     design_matrix=zeros(sum(numSubj)*numCond,sum(numSubj));
%     exchange_block = zeros(1,sum(numSubj)*numCond); 
%     row = 1;
%     diff_column = [];
%     diff_column1 = [];
%     diff_column2 = [];
%     for gg = 1:numGroups
%         for cc = 1:numCond
%             for ss = 1:ss
%                 % creates design matrix
%                 if gg == 1
%                     exchange_block(row) = ss ;
%                 else
%                     exchange_block(row) = ss + (sum(numSubj(1:gg-1)));
%                 end
%                 column = exchange_block(row);
%                 design_matrix(row,column) = 1;
%                 row = row + 1;
%                 % creats difference column in design_matrix
%                 if strcmp(network.difference,'Group')
%                     if gg == 1
%                         value = 1 ;
%                     elseif gg ==2
%                         value = -1 ;
%                     end
%                     diff_column = cat(1,diff_column,value);
%                 elseif strcmp(network.difference,'Condition')
%                     if cc == 1 
%                         value = 1 ;
%                     elseif cc ==2
%                         value = -1 ;
%                     end
%                     diff_column = cat(1,diff_column,value);
%                 elseif strcmp(network.difference,'Condition in Group')
%                     if gg == 1 % G1C2-G1C1
%                         if cc == 1
%                             value = -1 ;
%                         elseif cc ==2
%                             value = 1 ;
%                         end
%                         diff_column1 = cat(1,diff_column1,value);
%                     else
%                         value = 0;
%                         diff_column1 = cat(1,diff_column1,value);
%                     end
%                     if gg == 2 % G2C2-G2C1
%                         if cc == 1
%                             value = -1 ;
%                         elseif cc ==2
%                             value = 1 ;
%                         end
%                         diff_column2 = cat(1,diff_column2,value);
%                     else
%                         value = 0;
%                         diff_column2 = cat(1,diff_column2,value);
%                     end
%                 end
%             end
%         end
%     end
%     if strcmp(network.difference,'Condition in Group')
%         diff_column = cat(2,diff_column1,diff_column2);
%         UI.contrast.ui=[zeros(1,size(design_matrix,2)),UI.contrast.ui];
%     else
%         if UI.contrast.ui(1) == 1 && UI.contrast.ui(2) == -1 % contrast specified in diff_column
%             UI.contrast.ui=[zeros(1,size(design_matrix,2)),1];
%         elseif UI.contrast.ui(1) == -1 && UI.contrast.ui(2) == 1% INVERSE contrast specified in diff_column
%             UI.contrast.ui=[zeros(1,size(design_matrix,2)),-1];
%         end
%     end
%     design_matrix = cat(2,design_matrix,diff_column);
% else
%     % non-paired groups
%     design_matrix=zeros(sum(numSubj),numGroups);
%     for gg=1:numGroups
%         temp_numSubj=[0 numSubj];
%         rows = sum(temp_numSubj(1:gg))+1:sum(temp_numSubj(1:gg+1));
%         design_matrix(rows,gg) = 1;
%     end
% end
% save


%% run NBS in command line

if network.cmdline == 1
    global nbs;
    % exchange block (within subject comparison)
    if numCond > 1
        UI.exchange.ui = ['[',num2str(exchange_block),']'];
    else
        UI.exchange.ui='';
    end
    UI.perms.ui= num2str(UI.perms.ui);            
    UI.alpha.ui= num2str(UI.alpha.ui);   
    UI.contrast.ui= ['[',num2str(UI.contrast.ui),']'];
    UI.design.ui=[network.file_out,'/design_matrix.mat'];
    nbs_thresh = UI.thresh.min:UI.thresh.increment:UI.thresh.max;
    nbs_pval = nan(numFreqs,length(nbs_thresh));
    for ff= 1:numFreqs
        fprintf(['\nRuning command line NBS for ',num2str(network.filt_freqs(ff,1)),'-',num2str(network.filt_freqs(ff,2)),' Hz...\n\n']);
        for tt=1:length(nbs_thresh)
            UI.thresh.ui = num2str(nbs_thresh(tt)); 
            UI.matrices.ui=[network.file_out,'/NBS_',num2str(network.filt_freqs(ff,1)),'-',num2str(network.filt_freqs(ff,2)),'Hz.mat'];
            NBSrun(UI,[]);
            if ~isempty(nbs.NBS.pval)
                nbs_pval(ff,tt) = nbs.NBS.pval(1);
                if network.saveNBS == 1
                    save([network.file_out,'/NBS_',num2str(network.filt_freqs(ff,1)),'-',num2str(network.filt_freqs(ff,2)),'Hz_thresh-',num2str(nbs_thresh(tt)),'.mat'],'nbs','-mat');
                end
                close all
            end
        end
    end
    nbs_iterations = cat(1,nbs_pval,nbs_thresh);
    openvar('nbs_iterations');
    save([network.file_out,'/NBS_thr_iterations.mat'],'nbs_iterations','-mat');
    disp('Done running through NBS Thresholds.');
else
    global nbs;
    NBS;
end
