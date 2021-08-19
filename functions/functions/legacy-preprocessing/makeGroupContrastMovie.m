function [ ] = makeGroupContrastMovie( save_path, control, clinical, varargin )
%MAKEGROUPCONTRASTMOVIE Make differences connectivity matrix movie
% 
% makeGroupContrastMovie(save_path, control, clinical)
% makeGroupContrastMovie(save_path, control, clinical, 'ParameterName', _ParameterValue_)
%
% This function makes a video of the connectivity matrices of two
% contrasting groups, as they change over time. It plots the difference
% matrix alongside the group connectivity matrices, as well as an average
% amplitude plot to give a sense of where in the trial the video is at.
% 
% makeGroupContrastMovie automatically tries to encode the resulting video
% in H.264 format (in an mp4 file) using ffmpeg to save space. If this
% package is not available on the system, the video will remain as an
% uncompressed *.avi file.
%
% *Required Arguments:*
%   save_path - The location of the resulting video.
% 
%   control   - The connectivity matrix of the first / control group
%               This should have dimensions [sensors x sensors x samples]
% 
%   clinical  - The connectivity matrix of the second / clinical group
%               Must be exactly the same size as the control matrix
% 
% *Optional Arguments:*
%   wndX            - The width of the resulting video. 
%                       Default: 1000px
%                     The size of the connectivity matrices and amplitude
%                     plot automatically scale depending on the width of
%                     the video. The connectivity matrices will always fill
%                     the width of the video and will be square. The
%                     amplitude plot will take up the remaining vertical
%                     space.
%                     Note: The width must be divisible by 2
% 
%   wndY            - The height of the resulting video.
%                       Default: 720px
%                     Note: The height must be divisible by 2
% 
%   timeWindow      - The window of time (in samples) over which the video
%                     will animate.
%                       Default: the entire matrix
%                     This must be supplied as a 1x2 vector taking the form
%                     of [start, end].
% 
%   ampTimeWindow   - The window of time (in samples) that should be
%                     plotted on the amplitude plot.
%                       Default: the entire matrix
% 
%   controlAmp      - A custom amplitude time series to plot for the
%                     control group
%                       Default: the mean of group sensor pairs over time
% 
%   clinicalAmp     - A custom amplitude time series to plot for the
%                     clinical group
%                       Default: the mean of group sensor pairs over time
% 
%   title           - Title that is displayed at the top of the video for
%                     identification.
%                       Default: Group Contrast Movie
% 
%   cmax            - Color axis maximum value
%                       Default: 95th percentile of all pairs through time
%                     Supply cmax as an absolute connectivity value, not a
%                     percentile.
% 
%   cmin            - Color axis minimum value
%                       Default: 5th percentile of all pairs through time
% 
%   colormap        - Custom colormap for difference matrix
%                       Default: loads from difference_redblue_colormap.mat
%                     Supply colormap in the MATLAB colormap format. See
%                     'doc colormap'
% 
%   fps             - Frames per second in the resulting video
%                       Default: 15
% 
%   controlTitle    - Label name for control group. Appears in the legend
%                     and connectivity matrix plot title.
%                       Default: Controls
% 
%   clinicalTitle   - Label name for clinical group. Appears in the legend
%                     and connectivity matrix plot title.
%                       Default: Clinical
% 
%   indicatorLoc    - Location of event indicators (ie. stimulus onset) in
%                     number of samples.
%                       Default: no indicators
% 
%   indicatorColor  - Color of event indicators in MATLAB normalized RGB
%                     values.
%                       Default: [0.8, 0, 0]    (some shade of red)
%   
%   markerColor     - Color of 'current position in video' marker given in
%                     MATLAB normalized RGB values.
%                       Default: [0, 0.7, 0]    (some shade of green)
% 
%   ampXTickPos     - Custom X-axis tick positions for amplitude plot.
%                       Default: no change from MATLAB default
% 
%   ampXTickLabel   - Custom X-axis tick labels for amplitude plot.
%                       Default: no change from MATLAB default
%                     Note: Labels won't be set unless BOTH ampXTickPos and
%                     ampXTickLabel are both set.
% 
%   ampXLabel       - Custom X-axis title/label for amplitude plot
%                       Default: "Time (samples)"
%
% Simeon Wong
% 2013 July 5
%
% Edited: 2014 May 22

fprintf('%s - Generating group contrast movie...\n', datestr(now, 13));
tMovie = tic;

%% Process parameter arguments
p = inputParser;

addParameter(p, 'wndX', 1000);
addParameter(p, 'wndY', 800);
addParameter(p, 'timeWindow', []);
addParameter(p, 'ampTimeWindow', []);
addParameter(p, 'controlAmp', []);
addParameter(p, 'clinicalAmp', []);
addParameter(p, 'title', 'Group Contrast Movie');
addParameter(p, 'cmax', []);
addParameter(p, 'cmin', []);
addParameter(p, 'colormap', []);
addParameter(p, 'fps', 15);
addParameter(p, 'controlTitle', 'Controls');
addParameter(p, 'clinicalTitle', 'Clinical');
addParameter(p, 'indicatorLoc', []);
addParameter(p, 'indicatorColor', [0.8,0,0]);
addParameter(p, 'markerColor', [0,0.7,0]);
addParameter(p, 'ampXTickPos', []);
addParameter(p, 'ampXTickLabel', []);
addParameter(p, 'ampXLabel', 'Time (samples)');

parse(p, varargin{:});

%% Input matrix error checking
% Get size of input data matrices
sizeControl = size(control);
sizeClinical = size(clinical);

% Sanity check. Make sure input matrices match in size, and are square
if sizeControl(1) ~= sizeControl(2)
    error('Control input matrix is not a square connectivity matrix!');
end
if sizeClinical(1) ~= sizeClinical(2)
    error('Clinical input matrix is not a square connectivity matrix!');
end
if sizeControl(3) ~= sizeClinical(3)
    error('Connectivity matrices don''t have the same number of samples!');
end

num_frames = sizeControl(3);

%% Get difference matrix and amplitude waveforms

% Get difference matrix
difference = control - clinical;

% Calculate mean amplitudes (or use supplied parameters)
if isempty(p.Results.controlAmp)
    inveye = ~eye(sizeControl(1));                  % inverse identity matrix
    controlAmp = size(sizeControl(3),1);            % preallocate amplitude vector
    for kk=1:num_frames
        % Take the mean, of all the sensor pairs          v ensuring diagonal is 0
        controlAmp(kk) = mean(squareform(control(:,:,kk).*inveye, 'tovector'));
    end
else
    controlAmp = p.Results.controlAmp(:);
end

if isempty(p.Results.clinicalAmp)
    inveye = ~eye(sizeClinical(1));
    clinicalAmp = size(sizeClinical(3),1);
    for kk=1:num_frames
        clinicalAmp(kk) = mean(squareform(clinical(:,:,kk).*inveye, 'tovector'));
    end
else
    clinicalAmp = p.Results.clinicalAmp(:);
end

% Sanity check. Make sure amplitude waveforms have same number of samples
% as the number of frames.
if length(controlAmp) ~= num_frames || length(clinicalAmp) ~= num_frames
    error('Input amplitude waveforms do not have the same number of samples as the connectivity matrices!');
end

%% Process color mapping

% If color axis aren't explicitly specified, use 5th and 95th percentiles
if isempty(p.Results.cmax)
    cmax = prctile([control(:); clinical(:)], 95);
else
    cmax = p.Results.cmax;
end
if isempty(p.Results.cmin)
    cmin = prctile([control(:); clinical(:)], 5);
else
    cmin = p.Results.cmin;
end

diffcmax = prctile(difference(:), 95);
diffcmin = prctile(difference(:), 5);

% Make sure difference colour axes are centered about zero
diffcmax = max(abs(diffcmax), abs(diffcmin));
diffcmin = -1 * diffcmax;

% If custom colour axis for difference matrix not specified, use blue/red
if isempty(p.Results.colormap)
    load('difference_redblue_colormap.mat', 'cmap');
else
    cmap = p.Results.colormap;
end

%% Make movie

% temporary video save path
[fp_path, fp_name, ~] = fileparts(save_path);
temp_path = fullfile(fp_path, [fp_name, '.avi']);

% Initialize video writer for output
vid = VideoWriter(temp_path);
vid.Quality = 100;
vid.FrameRate = p.Results.fps;
open(vid);

figSize = [p.Results.wndX, p.Results.wndY];
% Ensure figure size is even, otherwise, ffmpeg whines about it
figSize(1) = figSize(1) + mod(figSize(1), 2);
figSize(2) = figSize(2) + mod(figSize(2), 2);

f = figure('Position', [10 50 figSize]);

% Positioning stuff
cbHeight = 60;          % Vertical height of colourbar
labelHeight = 25;       % Height of label area between matrices and amp plot

matBox = figSize(1) / 3;    % Size of bounding box for matrices
matSize = matBox - 40;      % Matrices plot size

ampHeight = figSize(2) - matBox - cbHeight - labelHeight - 80;      % Height of amp plot
ampWidth = figSize(1) - 80;                                         % Width of amp plot

% Store size and position, and convert to MATLAB normalized values
ampLoc = [60, matBox + cbHeight + labelHeight + 10, ampWidth, ampHeight] ./ [figSize, figSize];
matLoc = [20,          50+cbHeight, matSize, matSize;
          matBox+20,   50+cbHeight, matSize, matSize;
          matBox*2+20, 50+cbHeight, matSize, matSize;] ./ repmat(figSize,3,2);
      
% Set time window
if isempty(p.Results.timeWindow)
    timeWindow = [1, num_frames];
else
    timeWindow = p.Results.timeWindow;
end
if isempty(p.Results.ampTimeWindow)
    ampTimeWindow = [1, num_frames];
else
    ampTimeWindow = p.Results.ampTimeWindow;
end

for kk=timeWindow(1):timeWindow(2)
    figure(f);
    clf reset;
    
    % Plot amplitudes
    ax = subplot('Position', ampLoc);
    plot(controlAmp(ampTimeWindow(1):ampTimeWindow(2)));
    hold on
    plot(clinicalAmp(ampTimeWindow(1):ampTimeWindow(2)), 'r');
    
    xlabel(p.Results.ampXLabel);
    ylabel('Mean Connectivity');
    title(p.Results.title);
    
    legend(p.Results.controlTitle, p.Results.clinicalTitle);
    
    % Draw marker line
    v = axis;
    del = abs(v(4) - v(3)) * 0.05;
    ym = [v(3) + del, v(4) - del];
    line([kk,kk], ym, 'LineStyle', '--', 'Color', p.Results.markerColor);
    
    num_indicators = length(p.Results.indicatorLoc);
    for nn = 1:num_indicators
        line([1,1]*p.Results.indicatorLoc(nn)-ampTimeWindow(1)+1, ym, 'Color', p.Results.indicatorColor, 'LineStyle', ':');
    end
    
    % Set custom time axis tick positions
    if ~isempty(p.Results.ampXTickPos) && ~isempty(p.Results.ampXTickLabel)
        set(gca, 'XTick', p.Results.ampXTickPos - ampTimeWindow(1)+1, 'XTickLabel', p.Results.ampXTickLabel);
    end
    
    % Plot control matrix
    ax = subplot('Position', matLoc(1,:));      % create axis
    imagesc(control(:,:,kk));                   % plot matrix image
    caxis(ax, [cmin cmax]);                     % set consistent color axis
    colorbar('peer', ax, 'SouthOutside');       % show colorbar
    title(p.Results.controlTitle);              % set axis title
    axis square                                 % make plot square
    set(ax, 'Position', matLoc(1,:));           % move axis into position
    
    % Plot clinical matrix
    ax = subplot('Position', matLoc(2,:));
    imagesc(clinical(:,:,kk));
    caxis(ax, [cmin cmax]);                     % set consistent color axis
    colorbar('peer', ax, 'SouthOutside');
    title(p.Results.clinicalTitle);
    axis square
    set(ax, 'Position', matLoc(2,:));
    
    % Plot difference matrix
    ax = subplot('Position', matLoc(3,:));
    imagesc(difference(:,:,kk));
    caxis(ax, [diffcmin diffcmax]);
    colormap(ax, cmap);
    colorbar('peer', ax, 'SouthOutside');
    title('Difference');
    axis square
    set(ax, 'Position', matLoc(3,:));
    
    writeVideo(vid, getframe(f));               % save the frame to video
end

close(f);
close(vid);

fprintf('%s - Encoding video...\n', datestr(now, 13));
encodeH264(temp_path);

tMovie = toc(tMovie);
fprintf('%s - Movie making completed in %.2f seconds.\n', datestr(now, 13), tMovie);

end

