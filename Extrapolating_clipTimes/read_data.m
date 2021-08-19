function [movie_data] = read_data(data)

% This is a helper function for generate_markerTimes used to read
% behavioral (.PSY) files

%% read in PSY data
opts = delimitedTextImportOptions("NumVariables", 5);

%% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ":";

%% Specify column names and types
opts.VariableNames = ["VarName1", "VarName2", "STARTTueDec1316", "VarName4", "VarName5"];
opts.VariableTypes = ["string", "string", "string", "string", "string"];

%% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

%% Specify variable properties
opts = setvaropts(opts, ["VarName1", "VarName2", "STARTTueDec1316", "VarName4", "VarName5"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName1", "VarName2", "STARTTueDec1316", "VarName4", "VarName5"], "EmptyFieldRule", "auto");

%% Import the data
movie_data = readmatrix(data, opts);
movie_data = movie_data(:,4);

%% Clear temporary variables
clear opts
end