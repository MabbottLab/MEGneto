function labels = saveSensorsToFile(paths, subj_match, ss, output)

numView = 10;
cfg                         = [];
cfg.dataset                 = cellfun(@(x) [paths.rawdata '/' x], subj_match.ds, 'UniformOutput', false);
cfg.dataset                 = cfg.dataset{ss};
cfg.channel                 = {'MEG'};
sensor_level = ft_preprocessing(cfg);
cfg.channel                 = {'HLC0011','HLC0012','HLC0013', ...
    'HLC0021','HLC0022','HLC0023', ...
    'HLC0031','HLC0032','HLC0033'};
headpos = ft_preprocessing(cfg);

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
    'Position',[300 200 1500 700]);
% initial plot
for page_num = 1:numPages
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
            title([labels{range(index),1}],'fontweight','bold');
            set(gca,'XTickLabel','','XGrid','on');
        end
    end
    print(f,[output '/' num2str(page_num) 'page_channels.png'],'-dpng')
end
images = glob([output '/*page_channels.png']);
fullim = [];
for ii = 1:length(images)
    im = imread(images{ii});
    fullim = [fullim;im];
end
imwrite(fullim, [output '/' 'all_channels.png']) 
    