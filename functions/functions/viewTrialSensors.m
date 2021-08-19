function viewTrialSensors(p_strct, subj_location)

% set defaults
if nargin < 2 || isempty(subj_location)
    subj_location = 1; % default show first subject in p-structure
end

%%%%%%%%%%% import raw data
cfg                         = [];
cfg.dataset                 = p_strct.subj.subj_ds(subj_location);
cfg.channel                 = {'MEG'};
sensor_level = ft_preprocessing(cfg);
cfg.channel                 = {'HLC0011','HLC0012','HLC0013', ...
                               'HLC0021','HLC0022','HLC0023', ...
                               'HLC0031','HLC0032','HLC0033'};
headpos = ft_preprocessing(cfg);

% import new epoching data
load(p_strct.paths.trial_cfg(subj_location));

    
% prep work
fs = headpos.fsample;
numTrial = ceil(p_strct.tot_time / p_strct.epoch_period);

% remove part of data with 0mm HM (not possible... turned off MEG but still recording)
HM = headpos.trial{1,1}(1,:);
a = find(HM==0,1,'first');
if isempty(a)
    time = (1:sensor_level.sampleinfo(2))/fs; % sec
else
    sensor_level.trial{1,1}(:,a:end)=[];
    new_HM = headpos.trial{1,1}(:,1:a-1);
    time = (1:a-1)/fs;
end
data=sensor_level.trial{1,1};

% determine HM (head motion)
X = nan(length(time),9);    % [samples] x [grad channels]
for j = 1:9,
    x = 1000*new_HM(j,:);% in mm
    x = x-median(x); % median is almost equivalent to measured position relative to dewar (reference)
    X(:,j) = x;
end
y = max(abs(X),[],2);


%%%%%%%%%%%%% Make full view figure

% plot lines and text
vert_lines =[cfg.trl(:,1);cfg.trl(end,2)]./fs;
trl_text={};
for tt=1:numTrial
    trl_text=cat(1,trl_text,['Trial ',num2str(tt)]);
end

f1=figure('Visible','off','Units','pixels','Name','View Full Timeseries',...
    'Position',[300 200 1500 900]); 

% raw data plot
sp1=subplot(2,1,1);
plot(time,data);
title([p_strct.subj.ID{subj_location},' Raw MEG Data'],'fontweight','bold');
text((cfg.trl(:,1)./fs)+(0.12*p_strct.epoch_period),(max(data(:))+(std(data(:))/2)).*ones(numTrial,1),trl_text,'Parent',sp1); 
xlabel('Time (s)');
ylabel('Amplitude (V)');
axis([0 max(time) min(data(:))-(std(data(:))/2) max(data(:))+std(data(:))]);
vline(vert_lines,'r--');
set(gca,'XGrid','off'); 

% visualization of what HM and what time segments were created
sp2=subplot(2,1,2);
plot(time,y,'k');
title([p_strct.subj.ID{subj_location},' Head Motion with selected Trials'],'fontweight','bold');
text((cfg.trl(:,1)./fs)+(0.12*p_strct.epoch_period),(max(y(:))+(std(y(:))/2)).*ones(numTrial,1),trl_text,'Parent',sp2); 
ylabel('Head Motion (mm)');
xlabel('Time (s)');
axis([0 max(time) min(y(:))-(std(y(:))/2) max(y(:))+std(y(:))]);
vline(vert_lines,'r--');
set(gca,'XGrid','off'); 

hold on;
for i = 1:numTrial, % show segments made in green
    seg = cfg.trl(i,1):cfg.trl(i,2);
    plot(seg/fs,y(seg),'g');
end

set(f1,'Visible','on');

%%%%%%%%%% make GUI

f2=figure('Visible','off','Units','pixels','Name','View Each Trial Seperately',...
    'Position',[400 50 1500 900]); 

% initial plot
page_num=1;
range=cfg.trl(page_num,1):cfg.trl(page_num,2);
% raw data plot
subplot(2,1,1);
plot(time(range),data(:,range));
title([p_strct.subj.ID{subj_location},' Raw MEG Data: Trial ',num2str(page_num)],'fontweight','bold');
xlabel('Time (s)');
ylabel('Amplitude (V)');
axis tight;
set(gca,'XGrid','on'); 
% visualization of what HM and what time segments were created
subplot(2,1,2);
plot(time(range),y(range),'k');
title([p_strct.subj.ID{subj_location},' Head Motion: Trial ',num2str(page_num)],'fontweight','bold');
ylabel('Head Motion (mm)');
xlabel('Time (s)');
axis tight;
set(gca,'XGrid','on'); 

% buttons
btn_fwd=uicontrol('Style','pushbutton','String','>>>',...
    'Position', [400 20 120 20],...
    'Callback', @gen_trialplots_fwd);
btn_bwd=uicontrol('Style','pushbutton','String','<<<',...
    'Position', [200 20 120 20],...
    'Callback', @gen_trialplots_bwd);
% drop down menu
menu=cellstr((num2str([1:numTrial]','%d')));
popup=uicontrol('Style','popup','String',menu,...
    'Position',[1300 20 50 20],'Value',1,...
    'Callback', @jump_to_trialplot);
txt=uicontrol('Style','text',...
    'Position',[1220 20 80 20],...
    'String','Trial No.:');
set(f2,'Visible','on');
 
%%%%%%%%%%%%%%%%%%%%% go forward
function gen_trialplots_fwd (source,callbackdata)
    if page_num == numTrial % reached end - skip to begining
        page_num=1;
        popup.Value=page_num;
        range=cfg.trl(page_num,1):cfg.trl(page_num,2);
        % raw data plot
        subplot(2,1,1,'replace');
        plot(time(range),data(:,range));
        title([p_strct.subj.ID{subj_location},' Raw MEG Data: Trial ',num2str(page_num)],'fontweight','bold');
        xlabel('Time (s)');
        ylabel('Amplitude (V)');
        axis tight;
        set(gca,'XGrid','on'); 
        % visualization of what HM and what time segments were created
        subplot(2,1,2);
        plot(time(range),y(range),'k');
        title([(p_strct.subj.ID{subj_location},' Head Motion: Trial ',num2str(page_num)],'fontweight','bold');
        ylabel('Head Motion (mm)');
        xlabel('Time (s)');
        axis tight;
        set(gca,'XGrid','on'); 
    else % keep going forward
        page_num=page_num+1;
        popup.Value=page_num;
        range=cfg.trl(page_num,1):cfg.trl(page_num,2);
        % raw data plot
        subplot(2,1,1,'replace');
        plot(time(range),data(:,range));
        title([p_strct.subj.ID{subj_location},' Raw MEG Data: Trial ',num2str(page_num)],'fontweight','bold');
        xlabel('Time (s)');
        ylabel('Amplitude (V)');
        axis tight;
        set(gca,'XGrid','on'); 
        % visualization of what HM and what time segments were created
        subplot(2,1,2);
        plot(time(range),y(range),'k');
        title([p_strct.subj.ID{subj_location},' Head Motion: Trial ',num2str(page_num)],'fontweight','bold');
        ylabel('Head Motion (mm)');
        xlabel('Time (s)');
        axis tight;
        set(gca,'XGrid','on'); 
    end
end

%%%%%%%%%%%%%%%% go backward
function gen_trialplots_bwd(source,cllbackdata)
    if page_num == 1 % reached begining - skip to end 
        page_num=numTrial; % (last page)
        popup.Value=page_num;
        range=cfg.trl(page_num,1):cfg.trl(page_num,2);
        % raw data plot
        subplot(2,1,1,'replace');
        plot(time(range),data(:,range));
        title([p_strct.subj.ID{subj_location},' Raw MEG Data: Trial ',num2str(page_num)],'fontweight','bold');
        xlabel('Time (s)');
        ylabel('Amplitude (V)');
        axis tight;
        set(gca,'XGrid','on'); 
        % visualization of what HM and what time segments were created
        subplot(2,1,2,'replace');
        plot(time(range),y(range),'k');
        title([p_strct.subj.ID{subj_location},' Head Motion: Trial ',num2str(page_num)],'fontweight','bold');
        ylabel('Head Motion (mm)');
        xlabel('Time (s)');
        axis tight;
        set(gca,'XGrid','on'); 
    else % keep going backwards
        page_num=page_num-1;
        popup.Value=page_num;
        range=cfg.trl(page_num,1):cfg.trl(page_num,2);
        % raw data plot
        subplot(2,1,1,'replace');
        plot(time(range),data(:,range));
        title([p_strct.subj.ID{subj_location},' Raw MEG Data: Trial ',num2str(page_num)],'fontweight','bold');
        xlabel('Time (s)');
        ylabel('Amplitude (V)');
        axis tight;
        set(gca,'XGrid','on'); 
        % visualization of what HM and what time segments were created
        subplot(2,1,2,'replace');
        plot(time(range),y(range),'k');
        title([p_strct.subj.ID{subj_location},' Head Motion: Trial ',num2str(page_num)],'fontweight','bold');
        ylabel('Head Motion (mm)');
        xlabel('Time (s)');
        axis tight;
        set(gca,'XGrid','on'); 
    end
end
%%%%%%%%%%% jump to certain page
function jump_to_trialplot (source,callbackdata)
    page_num=source.Value;
    range=cfg.trl(page_num,1):cfg.trl(page_num,2);
    % raw data plot
    subplot(2,1,1,'replace');
    plot(time(range),data(:,range));
    title([p_strct.subj.ID{subj_location},' Raw MEG Data: Trial ',num2str(page_num)],'fontweight','bold');
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
    axis tight;
    set(gca,'XGrid','on'); 
    % visualization of what HM and what time segments were created
    subplot(2,1,2,'replace');
    plot(time(range),y(range),'k');
    title([p_strct.subj.ID{subj_location},' Head Motion: Trial ',num2str(page_num)],'fontweight','bold');
    ylabel('Head Motion (mm)');
    xlabel('Time (s)');
    axis tight;
    set(gca,'XGrid','on'); 
end


end