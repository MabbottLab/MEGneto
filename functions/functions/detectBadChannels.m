function [badChannels,resMat] = detectBadChannels(cfg,task)

thr = cfg.bchthr;
sections = cfg.sections;
dataset = cfg.dataset;

%if strcmp(task, 'Rest')
%    sample_epochs = cfg.epochs;
%else
%    disp('Non-Rest');
%end
%import raw data
cfg.dataset                 = cfg.dataset; %
cfg.channel                 = {'MEG'};
sensor_level = ft_preprocessing(cfg);
cfg.channel                 = {'HLC0011','HLC0012','HLC0013', ...
                               'HLC0021','HLC0022','HLC0023', ...
                               'HLC0031','HLC0032','HLC0033'};
headpos = ft_preprocessing(cfg);

% remove part of data with 0mm HM (not possible... turned off MEG but still recording)
HM = headpos.trial{1,1}(1,:);
a = find(HM==0,1,'first');
% cut off 2 secs from end (since there is noise as the MEG is turned off)
a = a-sensor_level.fsample*2;
cfg = [];
cfg.dataset                 = dataset; %
cfg.channel                 = {'MEG'};
cfg.lpfilter = 'yes';
cfg.lpfreq = [50];
cfg.hpfilter = 'yes';
cfg.hpfreq = [1];
% cfg.bpfilter = 'yes'; 
% cfg.bpfreq = [59 61];
% cfg.bpfiltord = p.loadds.bpfiltord; 
sensor_level2 = ft_preprocessing(cfg);

if ~isempty(a)
    sensor_level2.trial{1,1}(:,a:end)=[]; % remove end of recording
%    if strcmp(task, 'Rest')
%        sensor_level2.trial{1,1}(logical(repmat(sample_epochs(1:a-1),[size(sensor_level2.trial{1,1},1) 1]))) = NaN ; % set bad samples (artifact) to NaN
%    end
end

% Break trial up into sections 
data = sensor_level2; 
[r,c] = size(data.trial{1,1});
endV = floor(c/sections)*sections;
data.trial{1,1}(:,endV+1:end)=[];
[r,c] = size(data.trial{1,1});
data.trial{1,1} = permute(reshape(data.trial{1,1},[r,c/sections,sections]),[1,3,2]);

for i = 1:sections 
    dataTemp = sensor_level;
    dataTemp.trial{1,1} = squeeze(data.trial{1,1}(:,i,:));
    resMat(i,:) = detectionBadChannels_Algorithm(dataTemp,thr);
end

if any(any(resMat')) == 0
    disp('No bad channels detected');
    badChannels = [];
else
    [row,col] = find(resMat ==1);
    [n,b] = histc(col,unique(col));
    multiple = find(n>1);
    index = find(ismember(b,multiple));
    temp = col(index);
    ind = unique(temp);
    if ~isempty(ind) % more than 1 window is noisy
        for i=1:length(ind)
            badChannels{i} = data.label{ind(i)};      
        end
        disp('Bad Channel...');
        disp(badChannels);
    else
        disp('No bad channels detected');
        badChannels = [];
    end
end

end






