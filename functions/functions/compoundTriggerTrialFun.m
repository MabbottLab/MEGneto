function trl = compoundTriggerTrialFun(cfg)
%Defines a trial x from the prestim value (trigger - some specified value)
%to the start of the subsequent trial plus some poststim value (trigger x+1
%+ poststim value)
% Supports compound epoching of task related data. See epochingHelper.m for
% additional information on usage.
%
% Anne Keller

%% Function requirements %%
% Example trial function available at http://fieldtrip.fcdonders.nl/tutorial/preprocessing
% cfg.dataset
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue

%Abnormal Requirements (relative to some other fieldtrip epoching
%functions)
% cfg.trialdef.eventtype
% cfg.trialdef.eventlist
% cfg.trialdef.prestim

 
%get the header information from the data file
hdr = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);
trl = [];
name = '';
flag = 1;
trlNames = {};
for n = 1:size(cfg.trialdef.eventlist,2)
    trlNames{n,1} = n;
    trlNames{n,2} = cfg.trialdef.eventlist{n};
end
 

for i=1:length(event) %for each event, do the following
%     event(i).type
    [r,~] = find(strcmp(event(i).type, cfg.trialdef.eventtype));
    if r == 1
        if flag == 1
            begSample = event(i).sample + cfg.trialdef.prestim*hdr.Fs;
            endSample = event(i).sample + cfg.trialdef.poststim*hdr.Fs;
            offset = cfg.trialdef.prestim*hdr.Fs;
            name = event(i).type;
            flag = 0;
        else
            if isempty(find(strcmp(name, trlNames),1))
               trlNames{end+1,2} = name;
               trlNames{end, 1} = size(trlNames,1);
            end
            [triggerNumber, ~] = find(strcmp(name, trlNames));
            trl(end+1, :) = [round([begSample endSample]) offset  triggerNumber];
            begSample = event(i).sample + cfg.trialdef.prestim*hdr.Fs;
            endSample = event(i).sample + cfg.trialdef.poststim*hdr.Fs;
            offset = cfg.trialdef.prestim*hdr.Fs;
            name = event(i).type;
            
        end
    elseif i == length(event)
        [triggerNumber, ~] = find(strcmp(name, trlNames));
        trl(end+1, :) = [round([begSample endSample]) offset  triggerNumber]; 
    elseif ~isempty(find(strcmp(event(i).type, cfg.trialdef.eventtype), 1))
        name = strcat(name, event(i).type);
    end


    
end %end for

end

