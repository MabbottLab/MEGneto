function [ participant ] = load_participants( paths, step , prefix )
%LOAD_PARTICIPANTS A wrapper for readtable that uses the paths structure and has a
%default useful for Megne2.
if nargin < 3
    prefix = 'subj_';
end

fid = fopen(paths.([prefix step])); % open the file
participant = textscan(fid, '%s', 'Delimiter', '\n', 'whitespace', ''); % load in text using newline characters as separators
participant = cell2table(rmmissing(participant{1})); % remove the empty newline and convert to a table
fclose(fid); % close the file

if isempty(participant)
    error(['No participants specified - please add ds files to ' paths.conf_dir '/' prefix step '.csv']);
end
end