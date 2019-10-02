function [Class, FileDataset] = opencls(Dataset, Original, Verbose)
  
  % Class = opencls(Dataset, Original, Verbose)
  %
  % Generates a data structure array from ClassFile.cls in Dataset.  Each
  % element of the array (each trial class) contains the fields: Name, Id,
  % Count, Trials.  Trials is a vector of trial numbers, indexed from 1 of
  % length Count. It is allowed to give directly a .cls file instead of a
  % Dataset.
  %
  % Original [default false]: If true, opens the original class file, in
  % case it was modified previously.
  %
  % Verbose [default true]: Output information to the command window.
  %
  % Marc Lalancette, The Hospital for Sick Children, Toronto, Canada.
  % 2014-03-25

  % Parse input arguments.
  if ~exist('Dataset', 'var')
    Dataset = pwd;
  end
  if ~exist('Original', 'var')
    Original = false;
  end
  if ~exist('Verbose', 'var')
    Verbose = true;
  end
  
  % Allow specifying the file directly.
  if strcmpi(Dataset(end-3:end), '.cls')
    ClassFile = Dataset;
    Dataset = fileparts(Dataset);
  else
    ClassFile = [Dataset, filesep, 'ClassFile.cls'];
  end
  if ~exist(Dataset, 'file') || ~exist(ClassFile, 'file') || Original
    % See if there was a backup.
    ClassFileOriginal = [ClassFile(1:end-4), '_Original.cls'];
    if exist(ClassFileOriginal, 'file')
      if ~copyfile(ClassFileOriginal, ClassFile)
        error('Copy failed: %s to %s', ClassFileOriginal, ClassFile);
      end
      if Verbose
        fprintf('ClassFile.cls restored from %s.\n', ClassFileOriginal);
      end
    else
      error('Class file not found in %s.', Dataset);
    end
  end
  
  % Open file for reading.
  if Verbose
    fprintf('Opening trial class file: %s\n', ClassFile);
  end
  fid = fopen(ClassFile, 'rt', 'ieee-be');
  if (fid == -1)
    error('Failed to open file %s.', ClassFile);
  end
  
  % Prepare string formats.  This makes reading simpler.
  FileHeaderFormat = [ ...
    'PATH OF DATASET:\n', ...
    '%s\n\n\n']; % Dataset name with full path (starting with / on Linux).
  FileHeaderFormat2 = [ ...
    'NUMBER OF CLASSES:\n', ...
    '%d', ... % Number of trial classes.
    '\n\n\n']; % This brings us to beginning of first CLASSGROUPID line.
  
  % Read file header.
  FileDataset = fscanf(fid, FileHeaderFormat);
  %   FileDataset = char(FileDataset'); % Because there are mixed types, it was saved as column of char numbers.
  nClasses = fscanf(fid, FileHeaderFormat2);
  
  
  % Read marker data.
  if Verbose
    fprintf('Found classes ');
  end
  % Note: no use preallocating the structure since we don't know
  % how many trials we have for each marker.
  %   for m = 1:nClasses
  c = 1;
  while ~feof(fid)
    Field = fgetl(fid);
    switch Field(1:end-1)
      case 'NAME'
        Line = fgetl(fid);
        Class(c).Name = sscanf(Line, '%s', 1); %#ok<*AGROW>
        %       case 'COMMENT'
        %         Line = fgetl(fid);
      case 'CLASSID'
        Line = fgetl(fid);
        Class(c).Id = sscanf(Line, '%f', 1);
      case 'NUMBER OF TRIALS'
        Line = fgetl(fid);
        Class(c).Count = sscanf(Line, '%f', 1);
      case 'LIST OF TRIALS'
        % Get data
        fgetl(fid); % 'TRIAL NUMBER \n'
        % Trials are indexed from 0 in file, but in CTF programs they are
        % indexed from 1, so add 1 to match what users expect, as well as
        % Matlab indexing.
        Class(c).Trials = 1 + fscanf(fid, '%f', inf); % fscanf fills in column order.
        % Using inf instead of the expected Count brings us past the empty
        % lines, at the start of the next CLASSGROUPID line.  
        %
        % However, if Count is 0, it (sometimes? depending on Matlab
        % version?) also grabs the "C" of that line so we'd need to seek
        % back a character. But when it doesn't grab the "C", seeking back
        % would bring us back on an empty line and break the sequence.  And
        % since we don't care about CLASSGROUPID anyway, ignore this.
        %         if Markers(m).Count == 0
        %           fseek(fid, -1, 'cof');
        %         end
        if size(Class(c).Trials, 1) ~= Class(c).Count
          warning('Number of trials for class %s doesn''t match expected count.', Class(c).Name);
          Class(c).Count = size(Class(c).Trials, 1);
        end
        if Verbose
          fprintf('%s (%1.0f), ', Class(c).Name, Class(c).Count);
        end
        c = c + 1;
      case []
        % Safeguard for empty line.  Shouldn't happen, but if it does,
        % don't skip and just move on to the next line.
      otherwise
        fgetl(fid); % Skip other fields.
    end
  end
  
  fclose(fid);
  if Verbose
    fprintf('\b\b.\n\n');
  end
  
  if length(Class) ~= nClasses
    error('Expected %d classes, found %d.', nClasses, length(Class));
  end
  
  
end

% CTF Matlab code from readCTFds.m
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%   Function readClassFile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function TrialClass=readClassFile(ClassFileName);
%
% %   Reads a CTF ClassFile and stores information in a structure.
% %   The class file allows a user to store a list of classsified trials in a data set.
% %   The ClassFile format is defined in document CTF MEG File Formats, PN900-0088.
% %   This format is rigid and readClassFile assumes that the ClassFile has the format
% %   current in October 2006.
%
% %   Inputs :
% %      ClassFileName : marker file including the full path and extension .mrk.
% %      trialList : List of trials to read.  Trial numbering : 1,2,...
% %                  If omitted or empty, read all markers.
%
% %  Output : Structure array marker.  Output trial numbering starts at 1.
% %           See CTF MEG File Formats, (PN900-0088) for the meaning of the structure
% %           fields.  A trial mat=y start before the t=0 point, so it is possible to have
% %           markers with time<0 (see ds.res4.preTrigPts).
%
% TrialClass=struct([]);
%
% if exist(ClassFileName)~=2
%   return     % File doesn't exist.
% end
%
% fid=fopen(ClassFileName,'r','ieee-be');
%
% for k=1:5;fgetl(fid);end  % Skip 5 lines (including path info)
% nClass=sscanf(fgetl(fid),'%d',1); %Read value and skip to the start of the next non-blank line.
% if nClass<=0
%   fprintf('readClassFile: File %s has %d classes.\n',nClass);
%   return
% end
%
% TrialClass=struct('ClassGroupId',[],'Name',char([]),...
%   'Comment',char([]),'Color',char([]),'Editable',char([]),'ClassId',[],'trial',[]);
%
% for k=1:nClass
%   %  Find the start of the next class identification
%   %  There is no need to check for end of file because the loop ends before an attempt
%   %  is made to read class nClass+1.
%   while ~strcmp('CLASSGROUPID:',fgetl(fid));end
%   ClassGroupId=sscanf(fgetl(fid),'%d',1);
%   fgetl(fid);
%   Name=deblank(fgetl(fid));
%   fgetl(fid);
%   Comment=deblank(fgetl(fid));
%   fgetl(fid);
%   Color=deblank(fgetl(fid));
%   fgetl(fid);
%   Editable=deblank(fgetl(fid));
%   fgetl(fid);
%   ClassId=sscanf(fgetl(fid),'%d',1);
%   fgetl(fid);
%   No_of_Trials=sscanf(fgetl(fid),'%d',1);
%   fgetl(fid);fgetl(fid);
%   if No_of_Trials>0
%     trial=reshape(fscanf(fid,'%d',No_of_Trials),1,No_of_Trials);
%   else
%     trial=[];
%   end
%   %  Adjust trial numbering so it starts at 1.
%   TrialClass(k)=struct('ClassGroupId',ClassGroupId,'Name',Name,...
%     'Comment',Comment,'Color',Color,'Editable',Editable,'ClassId',ClassId,...
%     'trial',trial+1);
% end
% fclose(fid);
% end