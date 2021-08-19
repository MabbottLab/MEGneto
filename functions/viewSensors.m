function [pages] = viewSensors(rawdata, subj_location, numView)
try
    has_subj_ind = isempty(subj_location);
catch
    warning('subj_location poorly defined (viewSensors)')
    subj_location = 1;
end
% set defaults
if nargin < 3 || isempty(numView)
    numView = 10;          % default 10 sensor plots per page
elseif nargin < 2 || isempty(subj_location)
    subj_location = 1; % default show first subject in p-structure
end

%%%%%%%%%%% import raw data
cfg                         = [];
cfg.dataset                 = rawdata{subj_location};
cfg.channel                 = {'MEG'};
sensor_level = ft_preprocessing(cfg);
cfg.channel                 = {'HLC0011','HLC0012','HLC0013', ...
                               'HLC0021','HLC0022','HLC0023', ...
                               'HLC0031','HLC0032','HLC0033'};
headpos = ft_preprocessing(cfg);


% remove part of data with 0mm HM (not possible... turned off MEG but still recording)
HM = headpos.trial{1,1}(1,:);
a = find(HM==0,1,'first');
if isempty(a)
    time = (1:sensor_level.sampleinfo(2))/sensor_level.fsample; % sec
else
    sensor_level.trial{1,1}(:,a:end)=[];
    time = (1:a-1)/sensor_level.fsample;
end

data=sensor_level.trial{1,1};
labels=sensor_level.label;


%%%%%%%%%%%%% Make GUI

numPages=ceil(151/numView);

% generate page indices
pages(1,:)=[1,numView];
if mod(151,numView)==0
    for pp=2:numPages
        pages(pp,1)=pages(pp-1,1)+numView;
        pages(pp,2)=pages(pp-1,2)+numView;
    end
else
    for pp=2:numPages-1
        pages(pp,1)=pages(pp-1,1)+numView;
        pages(pp,2)=pages(pp-1,2)+numView;
    end
    pages(numPages,1)=(151-mod(151,numView)+1);
    pages(numPages,2)=151;
end

f=figure('Visible','off','Units','pixels','Name','Raw MEG Data From All Sensors',...
    'Position',[300 200 1500 900]); 

% initial plot
page_num=1;
range=pages(page_num,1):pages(page_num,2);
for index=1:length(range)
    if index == numView
        subplot(numView,1,index); % show time scale only on last subplot
        plot(time,data(range(index),:));
        title(labels{range(index),1},'fontweight','bold');
        xlabel('Time (s)');
        set(gca,'XGrid','on'); 
    else
        subplot(numView,1,index);
        plot(time,data(range(index),:));
        title(labels{range(index),1},'fontweight','bold');
        set(gca,'XTickLabel','','XGrid','on'); 
    end
end

% buttons
btn_fwd=uicontrol('Style','pushbutton','String','>>>',...
    'Position', [400 20 120 20],...
    'Callback', @gen_sensorplots_fwd);
btn_bwd=uicontrol('Style','pushbutton','String','<<<',...
    'Position', [200 20 120 20],...
    'Callback', @gen_sensorplots_bwd);
% drop down menu
menu=cellstr((num2str([1:numPages]','%d')));
popup=uicontrol('Style','popup','String',menu,...
    'Position',[1300 20 50 20],'Value',1,...
    'Callback', @jump_to_sensorplot);
txt=uicontrol('Style','text',...
    'Position',[1220 20 80 20],...
    'String','Page No.:');
set(f,'Visible','on');

%%%%%%%%%%%%%%%%%%%%% go forward
function gen_sensorplots_fwd (source,callbackdata)
    if page_num == numPages % reached end - skip to begining
        page_num=1;
        popup.Value=page_num;
        range=pages(page_num,1):pages(page_num,2);
        for index= 1:length(range)
            if index == numView % show time scale only on last subplot
                subplot(numView,1,index,'replace') ;
                plot(time,data(range(index),:));
                title(labels{range(index),1},'fontweight','bold');
                xlabel('Time (s)');
                set(gca,'XGrid','on'); 
            else
                subplot(numView,1,index,'replace') ;
                plot(time,data(range(index),:));
                title(labels{range(index),1},'fontweight','bold');
                set(gca,'XTickLabel','','XGrid','on'); 
            end
        end
    else % keep going forward
        page_num=page_num+1;
        popup.Value=page_num;
        range=pages(page_num,1):pages(page_num,2);
        for index= 1:length(range)
            if index == numView
                subplot(numView,1,index,'replace') % show time scale only on last subplot
                plot(time,data(range(index),:));
                title(labels{range(index),1},'fontweight','bold');
                xlabel('Time (s)');
                set(gca,'XGrid','on'); 
            else
                subplot(numView,1,index,'replace') 
                plot(time,data(range(index),:));
                title(labels{range(index),1},'fontweight','bold');
                set(gca,'XTickLabel','','XGrid','on'); 
            end
        end
        if length(range) ~= numView % last page (not enough plots) clear remaining plots
            for index = (length(range)+1):numView
                if index == numView
                    subplot(numView,1,index,'replace'); % show time scale only on last subplot
                    plot(time,data(1,:),'w'); % dummy plot to show time scale and keep subplot scale consistent
                    xlabel('Time (s)');
                    set(gca,'XGrid','on','Box','off','YTickLabel','','YTick',[]);
                else
                    subplot(numView,1,index,'replace');
                    set(gca,'XTickLabel','','XGrid','on','YTickLabel','','YTick',[]);
                end 
            end
        end
    end
end

%%%%%%%%%%%%%%%% go backward
function gen_sensorplots_bwd(source,cllbackdata)
    if page_num == 1 % reached begining - skip to end 
        page_num=numPages; % (last page)
        popup.Value=page_num;
        range=pages(page_num,1):pages(page_num,2);
        for index = 1:length(range)
           if index == numView
                subplot(numView,1,index,'replace') % show time scale only on last subplot
                plot(time,data(range(index),:));
                title(labels{range(index),1},'fontweight','bold');
                xlabel('Time (s)');
                set(gca,'XGrid','on'); 
            else
                subplot(numView,1,index,'replace') 
                plot(time,data(range(index),:));
                title(labels{range(index),1},'fontweight','bold');
                set(gca,'XTickLabel','','XGrid','on'); 
           end
        end
        if length(range) ~= numView % last page (not enough plots) clear remaining plots
            for index = (length(range)+1):numView
                 if index == numView
                    subplot(numView,1,index,'replace'); % show time scale only on last subplot
                    plot(time,data(1,:),'w'); % dummy plot to show time scale and keep subplot scale consistent
                    xlabel('Time (s)');
                    set(gca,'XGrid','on','Box','off','YTickLabel','','YTick',[]); 
                else
                    subplot(numView,1,index,'replace');
                    set(gca,'XTickLabel','','XGrid','on','YTickLabel','','YTick',[]);
                end 
            end
        end
    else % keep going backwards
        page_num=page_num-1;
        popup.Value=page_num;
        range=pages(page_num,1):pages(page_num,2);
        for index = 1:length(range)
            if index == numView
                subplot(numView,1,index,'replace') % show time scale only on last subplot
                plot(time,data(range(index),:));
                title(labels{range(index),1},'fontweight','bold');
                xlabel('Time (s)');
                set(gca,'XGrid','on'); 
            else
                subplot(numView,1,index,'replace') 
                plot(time,data(range(index),:));
                title(labels{range(index),1},'fontweight','bold');
                set(gca,'XTickLabel','','XGrid','on'); 
           end
        end
    end
end
%%%%%%%%%%% jump to certain page
function jump_to_sensorplot (source,callbackdata)
    page_num=source.Value;
    range=pages(page_num,1):pages(page_num,2);
    for index= 1:length(range)
        if index == numView
            subplot(numView,1,index,'replace') % show time scale only on last subplot
            plot(time,data(range(index),:));
            title(labels{range(index),1},'fontweight','bold');
            xlabel('Time (s)');
            set(gca,'XGrid','on'); 
        else
            subplot(numView,1,index,'replace') 
            plot(time,data(range(index),:));
            title(labels{range(index),1},'fontweight','bold');
            set(gca,'XTickLabel','','XGrid','on'); 
        end
    end
    if length(range) ~= numView % last page (not enough plots) clear remaining plots
        for index = (length(range)+1):numView
            if index == numView
                subplot(numView,1,index,'replace'); % show time scale only on last subplot
                plot(time,data(1,:),'w'); % dummy plot to show time scale and keep subplot scale consistent
                xlabel('Time (s)');
                set(gca,'XGrid','on','Box','off','YTickLabel','','YTick',[]);
            else
                subplot(numView,1,index,'replace');
                set(gca,'XTickLabel','','XGrid','on','YTickLabel','','YTick',[]);
            end 
        end
    end
end


end