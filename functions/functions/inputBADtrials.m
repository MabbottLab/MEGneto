function [goodTrials,badTrials,p]=inputBADtrials(p_strct,subj_location)

% prep work
numTrial = ceil(p_strct.epoching(subj_location,3) / p_strct.epoching(subj_location,2));

%%%%%% GUI
f=figure('Visible','off','Units','pixels','Name','Bad Trial Input',...
    'Position',[50 200 250 400]); 

% Enter Channels
textIN=uicontrol('Style','edit','String','Enter BAD trial number...',...
    'Position',[20 360 210 30],'Min',1,'Max',1,...
    'Callback', @add_to_list);


% List Channel Inputed
txt=uicontrol('Style','text',...
    'Position',[20 320 125 20],...
    'String','Bad Trial List :','FontWeight','bold');
textOUT=uicontrol('Style','text',...
    'Position',[20 60 210 255],'BackgroundColor','w');

% Finalize
btn_done=uicontrol('Style','pushbutton','String','Done',...
    'Position', [40 20 80 30],...
    'Callback', @gen_list);
btn_clear=uicontrol('Style','pushbutton','String','Clear',...
    'Position', [130 20 80 30],...
    'Callback', @clear_list);

set(f,'Visible','on');

% initial values
bad_list={};

function add_to_list (source,callbackdata)
    trial_num = num2str(source.String);
    bad_list=cat(1,bad_list,trial_num);
    textOUT.String = bad_list;
    textIN.String='';
end

function gen_list (source,callbackdata)
    % FIX (necessary???) - need to remove enough Trials so that subject has the same amount of Trials as everyone in study
%     if numTrial-length(bad_list) > needed_numTrial
%         warning('Not enough Trials removed.');
%     elseif numTrial-length(bad_list) < needed_numTrial
%         warning('Too many Trials removed.');
%     else
        disp('Done removing BAD Trials.');
        close all;
%     end
end

function clear_list (source,callbackdata)
    if length(bad_list) <= 1
        bad_list={};
    else
        bad_list=bad_list(1:end-1);
    end
    textOUT.String=bad_list;
end

% return outputs only after figures are closed
waitfor(f);
badTrials=str2double(bad_list)';
good_bin=ones(1,numTrial);
good_bin(badTrials)=0;
tmp=1:numTrial;
goodTrials=tmp(logical(good_bin));

% import new epoching data
load(p_strct.paths.epoching_info(subj_location,'mat'),'vectordata_text','markerSamples','y','time','fs','idxMaxOverallHM','maxOverallHM');
load(p_strct.paths.trl_cfg(subj_location),'cfg');
% reject bad trials by removing trials from variables (used in preprocessing)
markerSamples=markerSamples(goodTrials);
vectordata_text=vectordata_text(goodTrials,:); vectordata_text(:,1)=1:size(vectordata_text,1);
cfg.trl=cfg.trl(goodTrials,:);
% save updated variables
save(p_strct.paths.epoching_info(subj_location,'mat'),'vectordata_text','markerSamples','y','time','fs','idxMaxOverallHM','maxOverallHM','-mat');
save(p_strct.paths.trl_cfg(subj_location),'-mat');

% edit and save text file
formatSpec = 'Trial %d : %3.2f - %3.2f seconds\n';
file = fopen(p_strct.paths.epoching_info(subj_location,'txt'),'w');
fprintf(file,formatSpec,vectordata_text');
fclose(file);

% fix p-structure tot_time
p_strct.epoching(subj_location,3) = length(goodTrials)*p_strct.epoching(subj_location,2);
p = p_strct;
save(p_strct.paths.p_strct,'p','-mat');

% update visualization of what HM and what time segments were created 
figure('Name','Head Motion','Color',[1 1 1]);
set(gcf, 'Units', 'centimeters','PaperPositionMode','auto');
set(gcf, 'Position', [3 3 12 7]);

plot(time,y,'k'); % most of these variables were saved in epoching_#mmHM.mat
title(p_strct.subj.ID{ss},'interpreter','none');
ylabel('Head motion, mm','FontSize',12,'FontWeight','bold');
xlabel('Time, s','FontSize',12,'FontWeight','bold');
set(gca, 'FontSize', 12,'FontWeight','bold','Box','off');

hold on;
for i = 1:length(markerSamples), % show segments made in green
    seg = cfg.trl(i,1):cfg.trl(i,2);
    plot(seg/fs,y(seg),'g');
end
text(idxMaxOverallHM/fs - (1/fs), max(y), ['| Max HM: ', num2str(maxOverallHM), 'mm']);
saveas(gcf, p_strct.paths.fig_headmotion(subj_location));

    
end
