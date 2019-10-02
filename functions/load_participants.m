function [ participant ] = load_participants( paths, step , prefix )
%LOAD_PARTICIPANTS A wrapper for readtable that uses the paths structure and has a
%default useful for Megne2.
if nargin < 3
    prefix = 'subj_';
end
participant = readtable(paths.([prefix step]), 'ReadVariableNames', false, 'Delimiter', ',');
if isempty(participant)
    error(['No participants specified - please add ds files to ' paths.conf_dir '/' prefix step '.csv']);
end
end