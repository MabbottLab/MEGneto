function [goodChannelList,badChannelList] = inputBADchannels(rawdata,subj_location, subj_match)

% goodChannelList contain all MEG and reference channels except the MEG channels present in badChannelList
hdr=ft_read_header(rawdata(subj_location));
sel=false(size(hdr.label));
labels = hdr.label(cell2mat(cellfun(@(x) contains(x,'M'),hdr.label,'UniformOutput',false)));

f=figure('Visible','off','Units','pixels','Name','Bad Channel/Sensor Input',...
    'Position',[50 200 540 400]); 

ds_name = subj_match.ds{subj_location}; %You were working on this! Fix text entry UI
message = 'Enter BAD channel names...';
% Enter Channels
uicontrol('Style','text','String',ds_name,'Position',...
    [20 360 500 30]);
uicontrol('Style','text','String',message,'Position',...
    [20 330 500 30]);
textIN=uicontrol('Style','edit',...
    'Position',[20 50 250 230],'Min',5,'Max',5,...
    'Callback', @add_to_list);


% List Channel Inputed
txt=uicontrol('Style','text',...
    'Position',[20 280 500 20],...
    'String','Bad Channel List :','FontWeight','bold',...
    'Min',1,'Max',length(labels));
textOUT=uicontrol('Style','listbox','String', labels,...
    'Position',[280 50 250 230],'BackgroundColor','w');

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
    chan_name = upper(source.String);
    bad_list=cat(1,bad_list,chan_name);
    textOUT.String=bad_list;
    textIN.String='';
end

function gen_list (source,callbackdata)
    for i=1:length(bad_list);
        sel(strcmp(hdr.label,bad_list{i})) = 1;
    end
    if length(bad_list)~=nnz(sel)
        warning('Channel name not found in list... double check spelling');
    else
        disp('Done making BAD channel list.');
        close all;
    end
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
badChannelList=bad_list;
goodChannelList=hdr.label(~sel);

end