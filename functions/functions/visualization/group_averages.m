function [groupadjmat,diff_groupadjmat,numGroups,numDiff,numSources,filt_freqs] = group_averages(visual)

numGroups = size(visual.group_paths,2);

% load p-structures
for i=1:numGroups
    % load p-structure
    load(visual.group_paths{i},'p');
    % %%% INPUTS %%% (update anonymous function's workspace if any edits to subject paths made)
    p.paths.connmat =           @(id,postfix) [p.paths.output_dir(id), p.subj.ID{id}, '_', postfix,'.mat'];
    visual.p{i} = p;
end

% load data
matlinks = {};
for gg = 1:numGroups
    fprintf(['\n---- Group ',num2str(gg),' ----\n']);
    numSubj(gg) = nnz(visual.p{gg}.subj.include); % actual number of subjects per group (excluding subj)
    tmp_numSubj(gg) = length(visual.p{gg}.subj.ID); % all subjects in group
    matlinks_index = 1; %initialize
    for ss = 1:tmp_numSubj(gg)
        if visual.p{gg}.subj.include(ss) % only use included subjects
            fprintf(['\t---- Subject ',num2str(matlinks_index),' ----\n']);
            fprintf(['\t\t',visual.p{gg}.paths.connmat(ss,visual.connmetric),'\n']);
            matlinks{gg}{matlinks_index} = matfile(visual.p{gg}.paths.connmat(ss,visual.connmetric));                
            matlinks_index = matlinks_index + 1;
        end
    end
end

% get dimensions
m = matlinks{1}{1};
filt_freqs=visual.p{1}.filt_freqs;
[numSources, ~, ~, numFreqs] = size(m.adjmat);

for gg = 1:numGroups
    numSubj = length(visual.p{gg}.subj.ID); % can differ per group
    fprintf('Making averages for %s...\n', visual.group_names{gg});
    for ss = 1:numSubj
        temp_subjadjmat=matlinks{1, gg}{1, ss}.adjmat;
        [~, ~, numTrials, ~] = size(temp_subjadjmat); % can differ per subject
        % set diagonal to NaN, so that the mean isn't biased.
        eyenan = ones(numSources);
        eyenan(eye(numSources)==1) = NaN;
        temp_subjadjmat=temp_subjadjmat.*repmat(eyenan,1,1,numTrials,numFreqs); 
        % Average Trials [nchan nchan ntrial nfreq] -> [nchan nchan nfreq] 
        temp_subjadjmat=squeeze(mean(temp_subjadjmat,3,'omitnan'));
        groupadjmat(:,:,:,ss,gg)=temp_subjadjmat;  % [nchan nchan nfreq nsubj ngroup]
    end
end
% Rearrange matrix [nchan nchan nsubj nfreq ngroup]
groupadjmat=permute(groupadjmat,[1,2,4,3,5]); 
% Average Subj [nchan nchan nsubj nfreq ngroup] -> [nchan nchan nfreq ngroup]
groupadjmat=squeeze(mean(groupadjmat,3,'omitnan'));
% set diagonal to 1 (fully correlated)
groupadjmat(logical(repmat(eye(numSources),1,1,numFreqs,numGroups))) = 1;

% Group differences
numDiff= numGroups-1;
if numDiff~=0
    disp('Making differences...');
    diff_groupadjmat=zeros(numSources,numSources,numFreqs,numDiff);
    for dd = 1:numDiff
        diff_groupadjmat(:,:,:,dd)=groupadjmat(:,:,:,dd+1)-groupadjmat(:,:,:,dd);
    end
else
    diff_groupadjmat = [];
end

fprintf('Done!\n\n');

end