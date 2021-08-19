% This script is used to generate the marker times for clip and clippet
% changes. 

% The pseudocode below provides an outline of how this was
% completed:
% For each participant
%   For each run of data
%       Use the marker file to extract unique marker times from the marker
%       file into an array. Extract the sample difference between adjacent
%       marker time samples and identify the differences that are greater
%       than 400 as this is the treshold to identify clip changes.
%
%       Extract run-dependent clip presentation order from the .psy files
%       and load clip times which containg frame number of each clip and
%       clippet change (cliptimes_main.mat). Then, rearrange the clip time to
%       follow the presented clip order.
%
%       Combine marker and clip information:
%           Use the starting time sample of each clip to assign each marker
%           to a frame number. This will occur in increments of 10 as 
%           marker times increase in increments of 200 
%           (and 200 samples = 10 frames). Then, map the 
%           clip/clippet frame changes by adding or subtracting samples 
%           from the closest light trigger, as necessary. This will ensure 
%           that marker times are generated for each frame change in 
%           increments of 1 rather than 10. Save a table 
%           containing these clip marker times, movie order and clippet 
%           marker times. 

% SEE ALSO: read_data.m

%% Load paths
project_path = '/home/dmatic/OIRM_ppt_dataInfo'; % path where outputs will be stored 
MEG_path = '/home/dmatic/OIRM_PILOT';% path to OIRM data
fieldtrip_path = '/mnt/hpc-megneto/fieldtrip'; % path to the FieldTrip toolbox/functions 
addpath(fieldtrip_path)
ft_defaults;

%% Load ppt data
% extract participant IDs
prior_dir = pwd();
cd(MEG_path); %navigate to directory
ds_files = struct2table(dir('*.ds'));
ds_files = ds_files.name;
cd(prior_dir);

%% Need to manually remove the following from participant/files table
% ppt 05 - no behavioural data
% ppt 15 - no behavioural data
% ppt 18 - no behavioural data

%% Generate list of PIDs
PIDs = cell(1,length(ds_files));
for i = 1:length(ds_files)
    PIDs{i} = ds_files{i}(1:23);
end

% generate list of unique ppts
PIDs = unique(PIDs).';

%% Create directors for storing ppt data
for i = 1:length(PIDs)
    mkdir([project_path '/' PIDs{i}(1:3)]);
end

mkdir([project_path '/group']);

%% Pre allocate table for saving data
clipMarkers_allPpts = table();
clipMarks_allPpts_output = 'clipMarkers_allPpts.mat';

%% Code to generate the markers! 
for i = 1:length(PIDs) % for each ppt, get their data!! Remember each ppt has 1-2 MEG data files.
    
    %% Initialize some outputs
    markerDifferences  = 'markerDifferences.mat';
    markersOverlayed   = 'markersFigure.fig';
    
    %% Extract participant and enter their MEG folder to grab behavioral data (.psy file)
    current_MEGdata = contains(ds_files, PIDs{i});
    current_MEGdata = ds_files(current_MEGdata);
    
    prior_dir = pwd();
    cd([MEG_path '/Behavioural/' current_MEGdata{1}(1:3)]); % go to this participant's Behavioural data folder
    psy_file = dir('*.psy');
    cd(prior_dir);
    
    %% Behavioral data extraction
    subj_beh = [MEG_path '/Behavioural/' current_MEGdata{1}(1:3) '/' psy_file.name]; % get the path to the movie order data
    
    %% Go through all MEG data for the participant (there may be 2 runs for 1 ppt)
    for j = 1:length(current_MEGdata) 
        disp(current_MEGdata{j})
        subj_meg = [MEG_path '/' current_MEGdata{j}]; % get the path to the MEG data
        
        %% Generate events and markers
        % read events and isolate marker samples
        events = ft_read_event(subj_meg);
        events_cell = struct2cell(events)';
        marker_times = unique(cell2mat(events_cell(:,2))); % "marker_times" is equivalent to term "marker samples"

        % read channels
        dat = ft_read_data(subj_meg);
        subdat = dat(192:199,:);
        
        %% Extract differences b/w adjacent time samples and save it
        time_diff = diff(marker_times);
        save([project_path '/' current_MEGdata{j}(1:3) '/' current_MEGdata{j}(1:end-3) '_' markerDifferences],'time_diff','-v7.3')

        %% Read behavioural data
        movie_data = read_data(subj_beh);
        
        %% isolate movie order using the .psy file
        movie_data = movie_data(contains(movie_data,'/home'));
        movie_order = [];
        for k = 1:length(movie_data)
            movie_order = [movie_order, str2num(extractBetween(movie_data(k),47,48))];
        end
        
        %% Load clip times file
        clip_times = load('/media/mabbotthpf/datasets/OIRM_PILOT/scripts-data/cliptimes_main.mat');
        clips = clip_times.clipTimes(movie_order); % drop unused clips and re-order them
        
        % clarify the movie order for all particiapnts, some have funky #
        % of movies so have to go on a case by case basis for those 
        if matches(current_MEGdata{j}(1:3), '16_') % only had one run of data
            clips = clips;
            movie_order_tosave = movie_order;
        elseif matches(current_MEGdata{j}(1:3), 'P01') % only had 8 clips of data
            clips(6:7) = [];
            movie_order(6:7) = [];
            movie_order_tosave = movie_order;
        elseif contains(subj_meg, 'Free2') || contains(subj_meg, 'FreeViewing2') % for run 2 of clips we save clips 6-10
            clips = clips(6:10);
            movie_order_tosave = movie_order(6:10);
        else % for run 1 of clips we save clips 1-5
            clips = clips(1:5);
            movie_order_tosave = movie_order(1:5);
        end
        
        % if the marker time difference is less than 400 then it is not the
        % start of a clip (remember that data consists of clips and
        % clippets)
        if (time_diff(1) < 400) 
            marks = [time_diff(1)];
            index = [1];
        else
            marks = [];
            index = [];
        end
        
        for k = 1:length(time_diff)
            if (time_diff(k) > 400) 
                marks = [marks, time_diff(k)];% marks = the value in time_diff (so markers greater than 400)
                index = [index, find(time_diff == time_diff(k)) + 1]; % row where the spike greater than 400 occurs
                
            end
        end
        
        % participant 11 and 12 Free2 have uncessary 'marks' --> manually
        % drop them
        if matches(current_MEGdata{j},'10_MEG092_20161110_FreeViewing2.ds')
            marks = marks(1:5);
            index = index(1:5);
        elseif matches(current_MEGdata{j},'11_MEG085_20161115_Free2.ds')
            marks = marks(9:end);
            index = index(9:end);
        elseif (matches(current_MEGdata{j},'12_MEG085_20161209_Free2.ds'))
            marks = marks(4:end);
            index = index(4:end);
        elseif (matches(current_MEGdata{j}, 'P01_MEG092_20161031_Free.ds'))
            marks(6:17) = [];
            index(6:17) = [];
        end
        
        %% Indicate frame #s based on marker times
        for clip_start = 1:length(index) % each new 'c' is the start of a new movie
            if (clip_start == 8) && (matches(current_MEGdata{j}, 'P01_MEG092_20161031_Free.ds'))
                last = find(marker_times == marker_times(end, 1));
            elseif (clip_start == 5) && ~(matches(current_MEGdata{j}, 'P01_MEG092_20161031_Free.ds'))
                last = find(marker_times == marker_times(end, 1));
            else
                last = index(clip_start + 1) - 1;
            end

            for marker_rows  = index(clip_start):last
                if marker_rows == index(clip_start)
                    marker_times(marker_rows, 2) = 1;
                else
                    marker_times(marker_rows, 2) = marker_times(marker_rows - 1, 2) + 10;
                end
            end
        end

        %% Extrapolate clip changes
        adjusted_clips = cell(1, length(clips)); % pre-allocate space to store clip/clippet info
        
        for clip_num = 1:length(clips)
            if (clip_num == 8) && (matches(current_MEGdata{j}, 'P01_MEG092_20161031_Free.ds'))
                movie = marker_times(index(clip_num):find(marker_times == marker_times(end, 1)), :);
            elseif (clip_num == 5) && ~(matches(current_MEGdata{j}, 'P01_MEG092_20161031_Free.ds'))
                movie = marker_times(index(clip_num):find(marker_times == marker_times(end, 1)), :);
            else
                movie = marker_times(index(clip_num):index(clip_num+1) - 1, :);
            end

            data = movie(:,2); % only look at the second column of movie

            current_clip = zeros(1, length(clips{clip_num}));
            for clippet = 1:length(clips{clip_num})
                value = clips{clip_num}(clippet); 
                [val,idx] = min(abs(data - value));
                closestVal = data(idx);

                % find index number in the data
                index_value = find(movie(:,2) == closestVal);
                num_samples = (closestVal - value)*20;
                if closestVal > value
                    marker_sample = movie(index_value, 1) - num_samples;
                elseif closestVal < value
                    marker_sample = movie(index_value, 1) + num_samples;
                elseif closestVal == value
                    marker_sample = movie(index_value, 1);
                end

                current_clip(clippet) = marker_sample;
            end

            adjusted_clips{clip_num} = current_clip;
        end
        
        %% Save clip markers for participants 
        if (length(marks) == 5) || (matches(current_MEGdata{j}, 'P01_MEG092_20161031_Free.ds'))
            clip_marks = {{marker_times(index)}, {movie_order_tosave}, {adjusted_clips}};
            temp_table = cell2table(clip_marks, 'RowNames', {current_MEGdata{j}}, 'VariableNames', {'Clip Marks', 'Movie Order', 'Adjusted Clippet Changes'}); 
        end
        
        clipMarkers_allPpts = [clipMarkers_allPpts;temp_table];
        
        %% Plot clippet changes on light signal and save
        
        x = (marker_times(1, 1)-30):(marker_times(end, 1)+600); % look at 10 frames
        x = x(x >= 0);
        f = figure;
        k = 1;
        for m = [1 5:8] % plot every ADC channel
            subplot(5,1,k)
            plot(x, subdat(m,x))
            for n = 1:length(adjusted_clips)
                for o = 1:length(adjusted_clips{n})
                    xline(adjusted_clips{n}(o), '--r'); % clippet and clip changes
                end
            end
            hold off
            k = k + 1;
        end
        savefig(f, [project_path '/' current_MEGdata{j}(1:3) '/' current_MEGdata{j}(1:end-3) '_' markersOverlayed])
        close(f)
        
    end
    save([project_path '/group/' clipMarks_allPpts_output], 'clipMarkers_allPpts', '-v7.3');
end


