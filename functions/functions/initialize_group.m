function [p]=initialize_group(p)

% remove '/' from end of base path
temp_base_path = strsplit(p.paths.base,'/');
if isempty(temp_base_path{end})
    p.paths.base = strjoin(temp_base_path(1:end-1),'/');
end

% Path to subjects directories
[LIST, ISDIR] = glob([p.paths.base,'/*']);
% exclude groupANALYSIS directory 
if nnz(~cellfun(@isempty,strfind(LIST,'groupANALYSIS'))) > 0 
    exclude_dir= find(~cellfun(@isempty,strfind(LIST,'groupANALYSIS')));
    ISDIR(exclude_dir) = 0;
    LIST2 = LIST(ISDIR);
else
    LIST2 = LIST(ISDIR);
end
p.subj.subj_dir = LIST2;

% Path to subject raw MEG
% [LIST, ISDIR] = glob([p.paths.meg,'/*_Rest.ds']); % FIX
[LIST, ISDIR] = glob([p.paths.meg,'/*_',p.task,'.ds']);
LIST(~ISDIR);
%remove last character in string (/)
LIST=cellfun(@(x) x(1:end-1),LIST,'UniformOutput',false);
p.subj.subj_ds = LIST;

% Path to subject MRI file
[LIST, ISDIR] = glob([p.paths.mri,'/' p.mriformat]);
LIST(~ISDIR);
p.subj.subj_mri =  LIST;

% Subject ID's to run (Exact name of subject folder)
if ~isempty(strfind(p.paths.base,'*')) % subjects are found in sub-directories
    [LIST, ISDIR] = glob(p.paths.base); 
    if nnz(~cellfun(@isempty,strfind(LIST,'groupANALYSIS'))) > 0 
        exclude_dir= find(~cellfun(@isempty,strfind(LIST,'groupANALYSIS')));
        ISDIR(exclude_dir) = 0;
        LIST2 = LIST(ISDIR);
    else
        LIST2 = LIST(ISDIR);
    end
    temp_paths_base = LIST2; % sub-directories
    sub_dir=1;
    for i=1:length(p.subj.subj_dir)
        temp{i} = strsplit(p.subj.subj_dir{i},temp_paths_base{sub_dir,1});
        if length(temp{i}) < 2 % try next subdirectory
            sub_dir = sub_dir+1;
            temp{i} = strsplit(p.subj.subj_dir{i},temp_paths_base{sub_dir,1});
        end
        temp2{i} = temp{i}(2);
        newtemp{i} = strsplit(char(temp2{i}),'/');
        newtemp2{i} = newtemp{1,i}{1,1};
    end
else % subjects are found in directory
    for i=1:length(p.subj.subj_dir)
        temp{i} = strsplit(p.subj.subj_dir{i},[p.paths.base]);
        temp2{i} = temp{i}(2);
        newtemp{i} = strsplit(char(temp2{i}),'/');
        newtemp2{i} = newtemp{1,i}{1,2};
    end
end
p.subj.ID = newtemp2' ;
% disp(newtemp2);
% check if ds and mri files are of the same ID
sub_dir=1;
for i=1:length(p.subj.ID)
    
   try
       if ~isempty(strfind(p.paths.base,'*')) % subjects are found in sub-directories
           
           IDlength = size(char(p.subj.ID{i}),2);
           
           temp_dir = strsplit(p.subj.subj_dir{i},temp_paths_base{sub_dir,1});
           if length(temp_dir) < 2 % try next subdirectory
                sub_dir = sub_dir+1;
                temp_dir = strsplit(p.subj.subj_dir{i},temp_paths_base{sub_dir,1});
           end
           temp_dir = strrep(temp_dir{2},'/','');

           temp_ds = strsplit(p.subj.subj_ds{i},temp_paths_base{sub_dir,1});
           temp_ds = strsplit(temp_ds{2},'/');
           temp_ds = temp_ds{end};

           temp_mri= strsplit(p.subj.subj_mri{i},temp_paths_base{sub_dir,1});
           temp_mri = strsplit(temp_mri{2},'/');
           temp_mri = temp_mri{end};
           
       else % subjects are found in directory
  
           IDlength = size(char(p.subj.ID{i}),2);

           temp_dir = strsplit(p.subj.subj_dir{i},[p.paths.base]);
           temp_dir = strrep(temp_dir{2},'/','');

           temp_ds = strsplit(p.subj.subj_ds{i},[p.paths.base]);
           temp_ds = strsplit(temp_ds{2},'/');
           temp_ds = temp_ds{end}; %4

           temp_mri= strsplit(p.subj.subj_mri{i},[p.paths.base]);
           temp_mri = strsplit(temp_mri{2},'/');
           temp_mri = temp_mri{end}; %4
       end
   catch  % needed if missing file for last subj in list
        return % continue to next part
   end
             
   if ~strncmpi(temp_dir,p.subj.ID{i},IDlength) ...
           || ~strncmpi(p.subj.ID{i},temp_mri,IDlength) ...
           || ~strncmpi(p.subj.ID{i},temp_ds,IDlength)
       fprintf('Subject ID : %s does not correspond to one of the following:\n', char(p.subj.ID{i}) )
       fprintf('DS : %s \n', temp_ds);
       fprintf('MRI : %s \n', temp_mri);
        % check if DS exists for every subject
        if length(p.subj.subj_dir)~=length(p.subj.subj_ds)
            error(['Missing DS: ', p.subj.ID{i}]);
        elseif length(p.subj.subj_dir)~=length(p.subj.subj_mri)
        % check if MRI exists for every subject
            error(['Missing MRI:' , p.subj.ID{i}]);
        else
            error('Check files');
        end
   end
end

% check if head coil info present
for ds=1:length(p.subj.subj_ds)
    [LIST, ISDIR] = glob([p.subj.subj_ds{ds},'/*.hc']);
    LIST(~ISDIR);
%     disp(LIST);
    if ~isempty(LIST)
    [~,num_lines]=system(['grep -c ".$" ',LIST{1}]);
    num_lines=str2double(num_lines);
    if num_lines < 10
        error(['No Head Coil Positions Recorded: Remove subject', p.subj.ID{ds}]);
    end
    else
        error('Check if there are files in your .ds folder')
    end
end

%%% OUTPUTS
if ~isempty(strfind(p.paths.base,'*')) % subjects are found in sub-directories
    temp_base = strsplit(p.paths.base,'*');
    p.paths.group_output_dir = [temp_base{1,1},'groupANALYSIS/']; 
else
    p.paths.group_output_dir = [p.paths.base,'/groupANALYSIS/']; 
end
p.paths.p_strct = [p.paths.group_output_dir,'pINFO.mat']; 
p.paths.output_dir = @(id) [p.subj.subj_dir{id},'ANALYSIS/'];
%disp(p.paths.output_dir);


% Create Analysis folder for each subject if it doesn't exist
for ss= 1:length(p.subj.ID)
%     disp('Making Dir');
%     disp(ss);
%     disp(p.paths.output_dir(ss))
    if ~exist(p.paths.output_dir((ss)),'dir')
        mkdir(p.paths.output_dir((ss)));
    end
end

% Create Group Analysis folder if it doesn't exist
if ~exist(p.paths.group_output_dir,'dir')
    mkdir(p.paths.group_output_dir);    
end

% initialize subj.include parameter (default every subject is included)
p.subj.include = logical(ones(length(p.subj.ID),1));



% save p-structure in pINFO.mat
if exist(p.paths.p_strct,'file')
    % ask for user input to update info
    prompt = 'P-structure found... Do you want to overwrite it? [Y/N]:';
    %prompt = 'P-structure found... Do you want to: (Y: ) or (N:) [Y/N]:';
    answer = upper(input(prompt,'s'));
    if isempty(answer) || answer == 'Y' 
        disp('Updating p-structure.');
        tmp_p=p; % change name of p-structure info so it's not overwritten
        load(p.paths.p_strct,'p'); % load old p-structure
        % EPOCHING: update subject specific fields in p-structure (if applicable)
        if isfield(p,'epoching')
            for ss=1:length(p.subj.ID)
                tmp_index = strmatch(p.subj.ID{ss},tmp_p.subj.ID);
                if ~isempty(tmp_index)
                    present(ss,1)=tmp_index; % where old info present in new info
                    tmp_p.epoching(present(ss,1),:)=p.epoching(ss,:);
                end
            end
%             present % double check sorting
            disp('Updating epoching information.');
            if ~isempty(find(tmp_p.epoching(:,1)==0))
                warning('There is missing epoching information!');
                no_info_subj=find(tmp_p.epoching(:,1)==0);
                fprintf('Missing info for subject:\t%s\n',tmp_p.subj.ID{no_info_subj});
                fprintf('Options:\n\t1) Get info from other p-structures\n\t2) Manually enter the info and save the p-structure\n\t3) Re-run the specific subjects\n\n');
                prompt = 'What option do you choose? [1/2/3]:  ';
                answer = input(prompt);
                if answer == 1
                    stop_checking = 0;
                    while stop_checking == 0
                        prompt = 'Enter the path to a p-structure: ';
                        answer = input(prompt,'s');
                        try
                            check_p = load(answer,'p');
                            for ss=1:length(no_info_subj)
                                tmp_index = strmatch(tmp_p.subj.ID{no_info_subj(ss)},check_p.p.subj.ID);
                                if isempty(tmp_index)
                                    fprintf('[N]\tNo info found for subject %s.\n',tmp_p.subj.ID{no_info_subj(ss)});
                                else
                                    fprintf('[Y]\tInfo found for subject %s.\n',tmp_p.subj.ID{no_info_subj(ss)})
                                    present_info(ss,1)=tmp_index;
                                    tmp_p.epoching(no_info_subj(ss),:)=check_p.p.epoching(present_info(ss,1),:);
                                end
                            end
                        catch
                            warning(['No p-structure found in ', answer]);
                        end
                        prompt = 'Check any other p-structures?  [Y/N]: ';
                        answer = upper(input(prompt,'s'));
                        if answer == 'N'
                            stop_checking = 1;
                        end
                    end
                end
            end
            p.epoching=tmp_p.epoching;
            % save epoching info (just in case)
            epoch_info = p.epoching;
            save(p.paths.subj_epochInfo,'epoch_info', '-mat', '-v7.3');
        end
        % BAD COMPONENTS: update subject specific fields in p-structure (if applicable)
        if isfield(p,'bad_comp')
            for ss=1:length(p.subj.ID)
                tmp_index = strmatch(p.subj.ID{ss},tmp_p.subj.ID);
                if ~isempty(tmp_index)
                    present(ss,1)=tmp_index; % where old info present in new info
                    tmp_p.bad_comp{present(ss,1),1}=p.bad_comp{ss,1};
                end
            end
%             present % double check sorting
            disp('Updating ICA bad component information.');
            if  nnz(cellfun(@isempty,tmp_p.bad_comp)) > 0
                warning('There is missing ICA bad component information!');
                no_info_subj=find(cellfun(@isempty,tmp_p.bad_comp)); 
                fprintf('Missing info for subject:\t%s\n',tmp_p.subj.ID{no_info_subj});
                fprintf('Options:\n\t1) Get info from other p-structures\n\t2) Manually enter the info and save the p-structure\n\t3) Re-run the specific subjects\n\n');
                prompt = 'What option do you choose? [1/2/3]:  ';
                answer = input(prompt);
                if answer == 1
                    stop_checking = 0;
                    while stop_checking == 0
                        prompt = 'Enter the path to a p-structure: ';
                        answer = input(prompt,'s');
                        try
                            check_p = load(answer,'p');
                            for ss=1:length(no_info_subj)
                                tmp_index = strmatch(tmp_p.subj.ID{no_info_subj(ss)},check_p.p.subj.ID);
                                if isempty(tmp_index)
                                    fprintf('[N]\tNo info found for subject %s.\n',tmp_p.subj.ID{no_info_subj(ss)});
                                else
                                    fprintf('[Y]\tInfo found for subject %s.\n',tmp_p.subj.ID{no_info_subj(ss)})
                                    present_info(ss,1)=tmp_index;
                                    tmp_p.bad_comp{no_info_subj(ss),1}=check_p.p.bad_comp{present_info(ss,1),1};
                                end
                            end
                        catch
                            warning(['No p-structure found in ', answer]);
                        end
                        prompt = 'Check any other p-structures?  [Y/N]: ';
                        answer = upper(input(prompt,'s'));
                        if answer == 'N'
                            stop_checking = 1;
                        end
                    end
                end
            end
            p.bad_comp=tmp_p.bad_comp;
            % save all bad components (just in case)
            all_bad_comp = p.bad_comp;
            save(p.paths.ICAcomp_cfg,'all_bad_comp', '-mat', '-v7.3');
        end
        % BAD CHANNELS: update subject specific fields in p-structure (if applicable)
        if isfield(p,'bad_chann')
            for ss=1:length(p.subj.ID)
                tmp_index = strmatch(p.subj.ID{ss},tmp_p.subj.ID);
                if ~isempty(tmp_index)
                    present(ss,1)=tmp_index; % where old info present in new info
                    tmp_p.bad_chann{present(ss,1),1}=p.bad_chann{ss,1};
                end
            end
%             present % double check sorting
            disp('Updating bad channel information.');
            if  nnz(cellfun(@isempty,tmp_p.bad_chann)) > 0
                warning('There is missing bad channel information!');
                no_info_subj=find(cellfun(@isempty,tmp_p.bad_chann)); 
                fprintf('Missing info for subject:\t%s\n',tmp_p.subj.ID{no_info_subj});
                fprintf('Options:\n\t1) Get info from other p-structures\n\t2) Manually enter the info and save the p-structure\n\t3) Re-run the specific subjects\n\n');
                prompt = 'What option do you choose? [1/2/3]:  ';
                answer = input(prompt);
                if answer == 1
                    stop_checking = 0;
                    while stop_checking == 0
                        prompt = 'Enter the path to a p-structure: ';
                        answer = input(prompt,'s');
                        try
                            check_p = load(answer,'p');
                            for ss=1:length(no_info_subj)
                                tmp_index = strmatch(tmp_p.subj.ID{no_info_subj(ss)},check_p.p.subj.ID);
                                if isempty(tmp_index)
                                    fprintf('[N]\tNo info found for subject %s.\n',tmp_p.subj.ID{no_info_subj(ss)});
                                else
                                    fprintf('[Y]\tInfo found for subject %s.\n',tmp_p.subj.ID{no_info_subj(ss)})
                                    present_info(ss,1)=tmp_index;
                                    tmp_p.bad_chann{no_info_subj(ss),1}=check_p.p.bad_chann{present_info(ss,1),1};
                                end
                            end
                        catch
                            warning(['No p-structure found in ', answer]);
                        end
                        prompt = 'Check any other p-structures?  [Y/N]: ';
                        answer = upper(input(prompt,'s'));
                        if answer == 'N'
                            stop_checking = 1;
                        end
                    end
                end
            end
            p.bad_chann=tmp_p.bad_chann;
            % save all bad channels (just in case)
            all_bad_chann = p.bad_chann;
            save(p.paths.group_rmBadChan,'all_bad_chann', '-mat', '-v7.3');
        end
        % INCLUDED SUBJECTS: update subject specific fields in p-structure
        % -all subjects included by default, except those that were already
        % excluded in original p-structure
        for ss=1:length(p.subj.ID)
            tmp_index = strmatch(p.subj.ID{ss},tmp_p.subj.ID);
            if ~isempty(tmp_index)
                present(ss,1)=tmp_index; % where old info present in new info
                tmp_p.subj.include(present(ss,1))=p.subj.include(ss);
            end
        end
%             present % double check sorting
        % UPDATE OLD P-STRUCTURE
        p.paths.output_dir = tmp_p.paths.output_dir;
        p.paths.group_output_dir = tmp_p.paths.group_output_dir;
        p.paths.p_strct = tmp_p.paths.p_strct;
        p.subj= tmp_p.subj;
        p.task= tmp_p.task;
        p.mriformat = tmp_p.mriformat ;
        p.paths.base = tmp_p.paths.base;
        p.paths.meg = tmp_p.paths.meg;
        p.paths.mri = tmp_p.paths.mri;
        disp('Saving...');
        save(p.paths.p_strct,'p','-mat','-v7');
    else 
        disp('Not updating p-structure.');
    end
else
    disp('Saving...');
    save(p.paths.p_strct,'p','-mat','-v7');
end

disp('Done initializing subject paths and verifying files.');


end
