function [  ] = makeConnectivityMovie( savePath, adjmat, varargin )
%MAKECONNECTIVITYMOVIE Make a movie of the connectivity matrix over time,
%                      along with mean connectivity plot.
%
% Simeon Wong
% 2014 June 10
%
%%% Optional Arguments %%%
% timeWindow            - The start and end samples of the animation.
%                         Allows the user to show a longer time series, but
%                         only animate interesting portions of it.
%                           Default: The entire time series
% 
% fps                   - Frames Per Second: The number of frames per
%                         second to play the video. The higher this number,
%                         the "faster" the animation appears to go.
%                           Default: 15 frames per second
% 
% samplesPerFrame       - The number of samples averaged together to
%                         produce a single frame of the video. 
%                           Default: 1 (a one-to-one samples-to-video  
%                                    frame video)
% 
% connTitle             - Plot title for the connectivity matrix
%                           Default: (empty)
% 
% connXLabel            - x-axis label for connectivity matrix
%                           Default: 'Source Number'
% 
% connYLabel            - y-axis label for connectivity matrix
%                           Default: 'Source Number'
% 
% tsData                - Data to plot in the timeseries plot in the video
%                         to give an idea of where in the task we are
%                           Default: Mean of connectivity matrix
% 
% tsTitle               - Plot title for time series plot
%                           Default: (empty)
% 
% tsXLabel              - x-axis label for time series plot
%                           Default: 'Time (s)'
% 
% tsYLabel              - y-axis label for time series plot
%                           Default: 'Mean Connectivity'
% 
% cursorColor           - Colour of vertical line indicating current
%                         position in the time series that is shown in the
%                         connectivity matrix. Follows MATLAB ColorSpec.
%                           Default: 'r' (red)
% 
% cursorLineStyle       - Line style of cursor line. Follows MATLAB
%                         LineSpec.
%                           Default: '-' (solid line)
% 
% indicatorLoc          - Vector of sample numbers where indicator lines
%                         are drawn in the time series plot. This is useful
%                         for indicating onset of stimulus, etc..
%                           Default: no indicators
% 
% indicatorColor        - The color of vertical indicator
%                         lines. Supplied in MATLAB ColorSpec format.
%                           Default: [0, 0.7, 0] dark green
% 
% indicatorLineStyle    - The line style of the vertical indicator line.
%                           Default: dashed line ('--')
% 
% position              - Position and size of the figure window in a
%                         vector in pixels. [left, bottom, width, height]
%                           Default: [100 100 960 475]
% 
% connPosition          - Position and size of connectivity matrix plot
%                         within the figure window in pixels.
%                         [left, bottom, width, height]
%                           Default: [40 40 415 415]
%                           


%% Input parameters
p = inputParser;

addOptional(p, 'timeWindow', []);
addOptional(p, 'fps', 15);
addOptional(p, 'samplesPerFrame', 1);

addOptional(p, 'connTitle', '');
addOptional(p, 'connXLabel', 'Source Number');
addOptional(p, 'connYLabel', 'Source Number');

addOptional(p, 'tsData', []);
addOptional(p, 'tsTitle', '');
addOptional(p, 'tsXLabel', 'Time (s)');
addOptional(p, 'tsYLabel', 'Mean Connectivity');
addOptional(p, 'tsXTick', []);
addOptional(p, 'tsXTickLabel', {});

addOptional(p, 'cursorColor', 'r');
addOptional(p, 'cursorLineStyle', '-');

addOptional(p, 'indicatorLoc', []);
addOptional(p, 'indicatorColor', [0, 0.7, 0]);
addOptional(p, 'indicatorLineStyle', '--');

addOptional(p, 'position', [100 100 960 476]);
addOptional(p, 'connPosition', [40 40 415 415]);
addOptional(p, 'tsPosition', [515 40 435 200]);

parse(p, varargin{:});

%% Prepare inputs

% Take the mean of the connectivity matrix for meanall.
if isempty(p.Results.tsData)
    num_samples = size(adjmat, 3);
    meanall = zeros(num_samples,1);
    for kk = 1:num_samples
        meanall(kk) = mean(squareform(adjmat(:,:,kk), 'tovector'));
    end
else
    meanall = p.Results.tsData;
end

% Time window (in samples) within which to animate
if isempty(p.Results.timeWindow)
    timeWindow = [1 length(meanall)];
else
    timeWindow = p.Results.timeWindow;
end

% Color Axis
cx = [min(min(min(adjmat(:,:,timeWindow(1):timeWindow(2))))) max(max(max(adjmat(:,:,timeWindow(1):timeWindow(2)))))];

% y-axis for mean connectivity timeseries
ym = [min(meanall) max(meanall)];

% number of indicator lines
num_indicators = length(p.Results.indicatorLoc);

% temporary video save path
[fp_path, fp_name, ~] = fileparts(savePath);
tempPath = fullfile(fp_path, [fp_name, '.avi']);

%% Generate Movie

% Initialize video
vid = VideoWriter(tempPath);
vid.Quality = 100;
vid.FrameRate = p.Results.fps;
open(vid);

f = figure;
% Set size of figure
set(f, 'Position', p.Results.position);

for t=timeWindow(1):p.Results.samplesPerFrame:timeWindow(2)
    % Reset and clear the figure
    clf reset
    
    % Connectivity Matrix
    % Plot connectivity matrix
    subplot(1, 2, 1);
    imagesc(mean(adjmat(:,:,t*(1:p.Results.samplesPerFrame)),3));
    caxis(cx);
    axis square;
    
    % Set location of connectivity matrix
    set(gca, 'Units', 'pixels', 'Position', p.Results.connPosition);
    
    % Connectivity Matrix Labels
    xlabel(p.Results.connXLabel);
    ylabel(p.Results.connYLabel);
    title (p.Results.connTitle);
    
    % Time Series
    % Plot mean connectiity time series
    subplot(1, 2, 2);
    plot(meanall);
    line([t t], ym, 'Color', p.Results.cursorColor, 'LineStyle', p.Results.cursorLineStyle);
    
    % Set location of time series plot
    set(gca, 'Units', 'pixels', 'Position', p.Results.tsPosition);
    xlabel(p.Results.tsXLabel);
    ylabel(p.Results.tsYLabel);
    title (p.Results.tsTitle);
    
    % If set, set x-axis markers and labels
    if ~isempty(p.Results.tsXTick) && ~isempty(p.Results.tsXTickLabel)
        set(gca, 'XTick', p.Results.tsXTick, 'XTickLabel', p.Results.tsXTickLabel);
    end
    
    % Draw indicator lines
    for kk = 1:num_indicators
        line(p.Results.indicatorLoc(kk)*[1,1], ym, 'Color', p.Results.indicatorColor, 'LineStyle', p.Results.indicatorLineStyle);
    end
    
    % Get Frame
    set(f, 'PaperPositionMode', 'auto', 'Color', [1 1 1])
    refresh(f);
    
    writeVideo(vid, getframe(f));
end
close(f);
close(vid);

%% Encode video
encodeH264(tempPath);

end

