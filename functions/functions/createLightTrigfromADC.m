function createLightTrigfromADC(p,rangeOFsubj)
%
% Version 1 Sonya Bells
% 2016 November 
% Based on Mark L scripts
%
% INPUT FILES:
% - raw .ds (MEG 'RAW')
%
% OUTPUT FILES:
% - adjusts markerFile

%%% OUTPUTS %%%
% p.paths.fig_headmotion = @(id) [p.paths.output_dir(id),'headmotion.png'];


% determine which subjects to process
if nargin < 2 || isempty(rangeOFsubj)
    rangeOFsubj = 1:length(p.subj.ID);
else
    if strcmp(rangeOFsubj,'all')
        rangeOFsubj = 1:length(p.subj.ID);
    end 
end

% addpath('/data/Lily/megpipeline-Beamformer');
MatlabAddPaths;

for ss = rangeOFsubj

   fprintf('\n\n==================================\nSUBJECT: %s\n', p.subj.ID{ss});
   Dataset = p.subj.subj_ds{ss}; %p.subj.ID{ss};

    % New photodiode is in UADC005, 6 or 7, 8 (if forgot to change from arrows
    % task).
    %% Recreate light triggers ("OfflineLightOn/Off") from raw photodiode signal - ideally use the raw signal
    % channel 5
    % %% Recreate light triggers ("OfflineLightOn/Off") from filtered photodiode signal
    % channel 6
    %% Recreate light triggers ("OfflineLightOn/Off") from raw photodiode signal - ideally use the raw signal
    FindOptions.UseZeroLevel = false; % true for filtered 
    %FindOptions.AverageStartThresh = [-0.72, 0.98]; % use only if automatic finding of "start" didn't work.
    if isfield('p.trialdef','AverageStartThresh')
        FindOptions.AverageStartThresh = p.trialdef.light.AverageStartThresh; %[-0.72, 0.98]; % use only if automatic finding of "start" didn't work.
    end
    %FindOptions.AllowMissing = true;
    Overwrite = true;
        
    try
        MarkSignalOnOffEvents(Dataset, 5, FindOptions, Overwrite); % unfiltered new 
    catch
        warning('Light trigger error: Average all above threshold');
        continue
    end
    %% Check delay between Presentation codes and light trigger - There is a delay after stim code gets sent, so if using stim code, need to add offset, but if using photodiode, then no need
    MarkerTiming({Dataset}, {p.trialdef.t0marker1, p.trialdef.t0marker2}, 'OfflineLightOn');

    %% Add a marker for correct trials
    Mrk = openmrk(Dataset);
    tempMrk =  p.trialdef.details.includeOnce{end};
    TrialValidationFunction = @(x) (x(1) || x(2)) && x(3); %Look for event1 or event2 that has event 3 and add a new marker
    Mrk = AddMarker(Mrk, p.trialdef.newTrig, {p.trialdef.t0marker1, p.trialdef.t0marker2,tempMrk{1}, 'OfflineLightOn'}, ... %add 'Hit' to 'OfflineLightOn' for all events found
    TrialValidationFunction, 4, {p.trialdef.t0marker1, p.trialdef.t0marker2},Overwrite);
% 
    savemrk(Mrk, Dataset);

end
