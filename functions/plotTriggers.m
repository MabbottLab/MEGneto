function [numt0marker] = plotTriggers( dataset, t0marker , varargin )
%PLOTTRIGGERS Plot all the triggers in the file over time

p = inputParser;
addParameter(p, 'showFigure', false);
addParameter(p, 'savePath', []);
addParameter(p, 'xUnits', 'samples', @(x) any(strcmp(x, {'samples', 'sec', 'msec'})));
parse(p, varargin{:});

%% Processing

% extract list of events from the MEG dataset
eventslist = ft_read_event(dataset);

% get list of possible event types
uniquetypes = unique({eventslist.type});

% assign index to events (instead of just string)
eventslistnum = cellfun(@(str) find(strcmp(uniquetypes, str),1), {eventslist.type});

numt0marker = length(find( cellfun( @(X)isequal(X,t0marker),{eventslist.type})));
% length(find( cellfun( @(X)isequal(X,'RightCorrect'),{eventslist.type})))

%% Plot

xpoints = [eventslist.sample];

if strcmp(p.Results.xUnits, 'sec')
  hdr = ft_read_header(dataset);
  xpoints = xpoints / hdr.Fs;
elseif strcmp(p.Results.xUnits, 'msec')
  hdr = ft_read_header(dataset);
  xpoints = xpoints / hdr.Fs * 1000;
end

colourslist = lines(length(uniquetypes));
f = figure;
f.Position(3:4) = [1440 500];
for kk = 1:length(uniquetypes)
  plot(xpoints(eventslistnum == kk), kk * ones(sum(eventslistnum == kk),1), '+', 'LineWidth', 2, 'Color', colourslist(kk,:));
  hold on
end

ax = gca;
ax.YTickLabel = uniquetypes;
ax.YTick = 1:length(uniquetypes);

xlabel(['Time (', p.Results.xUnits, ')']);
ylabel('Trigger');

[~,filename,~] = fileparts(dataset);
title(filename, 'Interpreter', 'none');

if ~isempty(p.Results.savePath)
  saveas(f, p.Results.savePath);
end

if ~p.Results.showFigure
  close(f);
end

end

