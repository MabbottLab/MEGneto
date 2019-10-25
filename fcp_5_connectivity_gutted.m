function [stat] = fcp_5_connectivity_gutted(paths, source, source_proj)

% Gutted for pipeline workflow testing
disp('Just the illusion of connectivity analysis here.')
stat = sum(source.*source_proj);
pause(2)
disp('Quasi-done.')

end