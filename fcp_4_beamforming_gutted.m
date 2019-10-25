function [source, source_proj] = fcp_4_beamforming_gutted(paths)

% Gutted for pipeline workflow testing
disp('Just the illusion of beamforming going on here.')
source = randi(100,100,100);
source_proj = reshape(source(randperm(100.^2)),[100 100]);

end