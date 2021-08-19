function [  ] = megButterflyPlot( data, varargin )
%MEGBUTTERFLYPLOT 

%% parameters
p = inputParser;
addParameter(p, 'showFigure', false, @islogical);

parse(p, varargin{:});

%% plotting
f = figure;

end

