% HeadMotionTool
% Version 1.2 2016-03-11
% Marc Lalancette, The Hospital for Sick Children, Toronto, Canada.
%
% Display head movement information for a CTF MEG dataset collected with
% continuous head-coil localization.  The initial position saved in the
% dataset and typically intepreted as the unique head location, can be
% corrected to be the geometric median location, calculated from all trials
% (not marked as 'bad').  Trials can also be classified as 'bad' based on
% specified distance, range of motion and coil localization fitting error
% thresholds.  Position correction and trial rejection can be repeated
% since the geometric median location calculation does not include 'bad'
% trials.
%
% The function can be used with input arguments to allow scripting.  Every
% input is optional and can later be modified through the GUI.  Inputs must
% be provided as 'ArgumentName'/ArgumentValue pairs.
%
% [Location, MaxDistance] = HeadMotionTool('ArgName1', ArgValue1, ...)
%
%
% This tool can also be used with Fieldtrip, however in that case the CTF
% dataset is not modified in any way.  A "ft_definetrial" configuration
% structure (cfg) must be provided as input argument and two additional
% outputs are returned with rejected trials removed, and with corrected
% sensor positions:
%
% [Location, MaxDistance, cfg, grad] = HeadMotionTool('Fieldtrip', cfg, ...)
%
%
%  *Input arguments* (Default values in brackets):
%
% Dataset []: Full path to dataset directory.
%
% UseInitial [false]: If true, gives information based on dataset head coil
% information (usually initial position), instead of the calculated
% geometric median location.
%
% RejectThreshold [5]: Highlight trials as bad ('BAD_HeadMotion_Distance')
% based on the coil distances from the chosen reference location (median or
% initial), in mm. Depends on the proportion of samples within the trial
% that are above this threshold and the selected TrialRejectProportion (see
% below).  The range of motion within a trial is also considered and if it
% exceeds this threshold for more than the selected proportion, the trial
% is marked as bad ('BAD_HeadMotion_Range'). E.g. with a threshold of 5 mm,
% movement within one trial from -4 to +4mm in one direction would give a
% range of 8mm and that trial would be highlighted.
%
% RejectProportion [10]: Trials containing this proportion (in %) or more
% of samples exceeding the thresholds will be marked as bad.
%
% RejectFitThreshold [10]: Trials where the coil localization fitting error
% (greatest among the 3 coils) is above this threshold are marked as bad
% ('BAD_HeadMotion_Fit').
%
% CorrectInitial [false]: Modify Dataset (or Fieldtrip output structure) to
% use the calculated geometric median location to replace the "initial"
% position.  This allows for more accurate source localization since often
% this position is interpreted as *the* head position.  For CTF datasets
% (not Fieldtrip), this is done with the changeHeadPos CTF command line
% program, only available on Linux platforms.  This will only be done if
% UseInitial is false.
%
% SavePictureFile ['']: Filename (with or without .png extension) for
% saving a screenshot of the figure.  This is required to be able to
% automatically reject trials from the command line with the RejectTrials
% option.
%
% RejectTrials [false]: 
%
% Renderer ['painters']: Matlab option for the figure.  'opengl' allows
% transparency and smooth lines, but can be sometimes buggy and fonts are
% not smooth, while 'painters' or 'zbuffer' pretty much always work but
% don't allow transparency or smooth lines.
%
% GUI [true]: Set to false to prevent the GUI from appearing and use the
% function as a command line tool, useful for scripting.  Without the GUI,
% only initial position correction is possible, no trial rejection based on
% movement thresholds.  Using the GUI increases user awareness of head
% movement, which is a primary goal of this program.
%
%
%  *Ouptuts* (only possible if providing a valid dataset from the command line):
%
% Location: Position in dewar coordinates in mm, of the 3 coils:
%   [Na_x, y, z, LE_x, y, z, RE_x, y, z]
% Depending on UseInitial, either the median or initial location.
%
% MaxDistance: Maximum coil distance from Location (out of the 3 coils)
% across all trials and time samples.

%  *Dependencies* (other than functions in private subdirectory):
% For saving screenshot with 'SavePictureFile' option.
% For CTF datasets only:
%  1. Unofficial CTF Matlab package from MISL, provided with SPM or
%  Fieldtrip, last updated 2012-04-16 (Linux specific bug fix);
%  2. opencls, savecls and BackupOriginal for opening and saving trial
%  classification files; openhc for restoring datasets.
% For Fieldtrip only (other than Fieldtrip itself):
%  ChangeCoordinates

% This was an output argument before, now uiwait is in this function, the
% figure is now modal:
% FigureHandle: Useful for example for scripting and waiting until the
% HeadMotionTool figure is closed before moving on to the next dataset,
% using uiwait(FigureHandle).

function varargout = HeadMotionTool(varargin)
  
  %HEADMOTIONTOOL M-file for HeadMotionTool.fig
  %      HEADMOTIONTOOL, by itself, creates a new HEADMOTIONTOOL or raises the existing
  %      singleton*.
  %
  %      H = HEADMOTIONTOOL returns the handle to a new HEADMOTIONTOOL or the handle to
  %      the existing singleton*.
  %
  %      HEADMOTIONTOOL('Property','Value',...) creates a new HEADMOTIONTOOL using the
  %      given property value pairs. Unrecognized properties are passed via
  %      varargin to HeadMotionTool_OpeningFcn.  This calling syntax produces a
  %      warning when there is an existing singleton*.
  %
  %      HEADMOTIONTOOL('CALLBACK') and HEADMOTIONTOOL('CALLBACK',hObject,...) call the
  %      local function named CALLBACK in HEADMOTIONTOOL.M with the given input
  %      arguments.
  %
  %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
  %      instance to run (singleton)".
  %
  % See also: GUIDE, GUIDATA, GUIHANDLES
  
  % Last Modified by GUIDE v2.5 10-Mar-2016 11:49:54
  
  % Begin initialization code - DO NOT EDIT
  gui_Singleton = 1;
  gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @HeadMotionTool_OpeningFcn, ...
    'gui_OutputFcn',  @HeadMotionTool_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
    'gui_Callback',   []);
  if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
  end
  
  % Need to know nargout in opening function.  But figure doesn't exist yet
  % (can't use appdata).  Must set "global" variable here (or make custom
  % gui_mainfcn).
  
  if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
  else
    gui_mainfcn(gui_State, varargin{:});
  end
  % End initialization code - DO NOT EDIT
end


% --- Executes just before HeadMotionTool is made visible.
function HeadMotionTool_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
  % This function has no output args, see OutputFcn.
  % hObject    handle to figure
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  % varargin   unrecognized PropertyName/PropertyValue pairs from the
  %            command line (see VARARGIN)
    
  % See if dependencies should be added to the Matlab path (if Fieldtrip or
  % SPM are not found for example).  When this tool is distributed
  % "standalone" all additional dependencies should be in the private
  % directory, otherwise, as part of the larger repository, all required
  % functions should already be in the path.
  Dir = fileparts(mfilename('fullpath'));
  if ~exist('readCTFds.m', 'file') && ...
      isdir([Dir, filesep, 'FieldtripExternalCTF']) % "standalone"
    addpath([Dir, filesep, 'FieldtripExternalCTF']);
  elseif ~exist('readCTFds.m', 'file') || ~exist('opencls.m', 'file') || ...
      ~exist('ChangeCoordinates.m', 'file')
    if isdir([Dir, filesep, '..', filesep, 'External', filesep, 'ctf']) % in repository
      fprintf('Using MatlabAddPaths in repository to add all required folders.\n');
      CurDir = cd;
      cd([Dir, filesep, '..']);
      MatlabAddPaths; % May give error if repository is not complete.
      cd(CurDir);
    else
      error('Missing dependencies.  Ensure you have the full repository or the standalone version.');
    end
  end
  
  % Check if GUI should be displayed.  If not, do NOT load previous
  % parameters.  That would lead to unpredictable command line behavior and
  % render defaults useless.
  handles.ShowGUI = true;
  for a = 1:2:(nargin - 3),
    switch lower(varargin{a})
      case 'gui'
        if ~isempty(varargin{a+1})
          handles.ShowGUI = varargin{a+1};
        end
    end
  end

  % First load previous parameters if parameters file found in current
  % directory with default name.  These might get overwritten by
  % command-line parameters later.
  Dir = cd;
  handles.UserDataFile = [Dir, filesep, 'HeadMotionTool.mat'];
  if handles.ShowGUI && exist(handles.UserDataFile, 'file')
    load(handles.UserDataFile);
    % Load user parameters in text boxes.
    for EditBox = (fieldnames(UserData.Edit))'
      set(handles.(EditBox{1}), 'String', UserData.Edit.(EditBox{1}));
    end
    
    if UserData.uipanel_Ref
      set(handles.uipanel_Ref, 'SelectedObject', handles.radio_Initial);
    else
      set(handles.uipanel_Ref, 'SelectedObject', handles.radio_Median);
    end
    
    % Get directory of last loaded dataset.
    handles.DatasetDir = UserData.DatasetDir;
    % Other options.
    handles.FontSize = UserData.FontSize;
    set(handles.Menu_Opengl, 'Checked', UserData.Menu_Opengl);
    if strcmpi(UserData.Menu_Opengl, 'on')
      set(hObject, 'Renderer', 'opengl');
    else
      set(hObject, 'Renderer', 'painters');
    end
  else
    % Default values.
    handles.DatasetDir = '';
    handles.FontSize = 8;
    
    set(handles.uipanel_Ref, 'SelectedObject', handles.radio_Median);
    set(handles.edit_Thresh, 'String', '5.0'); % mm
    set(handles.edit_InterCoilThresh, 'String', '0.5'); % mm
    set(handles.edit_Proportion, 'String', '10.0');
    set(handles.edit_FitThresh, 'String', '10.0');
    set(hObject, 'Renderer', 'painters'); % 'opengl', 'zbuffer', 'painters'
    set(handles.Menu_Opengl, 'Checked', 'off');
  end

  
  % ___________________________________
  % Constants and other global and default values.
  handles.NumFormat = '%1.1f';
  handles.NumFormatSmall = '%1.2f';
  handles.MinDist = 0.01; % Must not be smaller than NumFormatSmall can display.
  if str2double(num2str(handles.MinDist, handles.NumFormatSmall)) == 0
    error('MinDist must not be smaller than NumFormatSmall can display.');
  end
  handles.NumFormatInt = '%1.0f';
  handles.Phi = (1 + sqrt(5))/2;
  handles.DefaultFigurePosition = [520, 400, 600, 400];
  handles.DefaultAxisPosition = [40, 60, 520, 180];
  %   handles.DefaultLeftPanelPosition = [10, 280, 190, 120];
  %   handles.DefaultMidPanelPosition = [210, 250, 180, 150];
  %   handles.DefaultRightPanelPosition = [400, 250, 190, 150];
  %   handles.DefaultOverlayPanelPosition = [210, 250, 380, 150];
  %   handles.DefaultButtonPosition = [10, 250, 120, 23];
  handles.PanelPositions = handles.DefaultFigurePosition(4) - [250, 280]; % From figure top.
  handles.AxisTop = handles.DefaultFigurePosition(4) - ...
    handles.DefaultAxisPosition(2) - handles.DefaultAxisPosition(4); % From figure top.
  
  handles.TextOptions = {'FontSize', handles.FontSize, ... % 'FontName', 'Helvetica',
    'Units', 'pixels', 'Margin', 1, 'Clipping', 'off', ...
    'HorizontalAlignment', 'center'};
  handles.TextPos = 37;
  
  % For screenshot, need to offset by window frame. (trial and error
  % depending on platform...)
  if ispc % Windows
    handles.WindowFrameSize = 8; % Pixels
  else % Linux (mac not tested)
    handles.WindowFrameSize = [2, 4]; % Pixels
  end
  
  % Backup of Proportion value.
  handles.Proportion = 10;

  % Defaults.
  RejectTrials = false;
  CorrectInitial = false;
  handles.Fieldtrip = false;
  handles.SavePictureFile = '';
  
  
% ___________________________________
  % Get and test input arguments.  Some default values are filled above,
  % (some only when there was no UserData), others are saved in the figure
  % itself.
  if isodd(nargin - 3)
    error('HeadMotionTool:OddInput', 'Odd number of inputs.  They must be provided in ''Name''/Value pairs.');
  end
  for a = 1:2:(nargin - 3),
    switch lower(varargin{a})
      case 'dataset'
        Dataset = varargin{a+1};
      case 'usefids'
        UseInitial = varargin{a+1};
      case 'correctinitial'
        if ~isempty(varargin{a+1})
          CorrectInitial = varargin{a+1};
        end
      case 'rejecttrials'
        if ~isempty(varargin{a+1})
          RejectTrials = varargin{a+1};
        end
      case 'rejectthreshold'
        RejectThreshold = varargin{a+1};
      case 'rejectproportion'
        RejectProportion = varargin{a+1};
      case 'rejectfitthreshold'
        RejectFitThreshold = varargin{a+1};
      case 'coilthreshold'
        InterCoilThreshold = varargin{a+1};
      case 'savepicturefile'
        handles.SavePictureFile = varargin{a+1};
      case 'renderer'
        Renderer = varargin{a+1};
      case 'gui' % Now checked earlier.
        %handles.ShowGUI = varargin{a+1};
      case 'fieldtrip'
        handles.Fieldtrip = true;
        handles.ft.cfg = varargin{a+1};
        handles.ft.grad = struct(); % If can't load data for some reason, to avoid output function error.
      otherwise
        error('HeadMotionTool:UnknownInput', 'Unrecognized input: %s.', varargin{a});
    end
  end
  
  % Check if correction is possible.
  if handles.Fieldtrip
    handles.CorrectAvailable = true;
  elseif isunix 
    [Status, Output] = system('changeHeadPos --help'); %#ok<NASGU> % Output required to suppress output in Matlab command window.
    if ~Status
      handles.CorrectAvailable = true;
    else
      handles.CorrectAvailable = false;
    end
  else
    handles.CorrectAvailable = false;
  end
  
  % Setting a radio button programmatically (through button or buttongroup)
  % respects exclusivity, but it doesn't trigger the group's selection
  % change function.
  if ~exist('UseInitial', 'var') || isempty(UseInitial)
  elseif UseInitial
    set(handles.uipanel_Ref, 'SelectedObject', handles.radio_Initial);
  end
  
  if ~exist('RejectThreshold', 'var') || isempty(RejectThreshold)
  elseif RejectThreshold < 0
    error('RejectThreshold must be positive.');
  else
    set(handles.edit_Thresh, 'String', num2str(RejectThreshold, handles.NumFormat));
  end
  if ~exist('InterCoilThreshold', 'var') || isempty(InterCoilThreshold)
  elseif InterCoilThreshold < 0
    error('InterCoilThreshold must be positive.');
  else
    set(handles.edit_InterCoilThresh, 'String', num2str(InterCoilThreshold, handles.NumFormat));
  end
  if ~exist('RejectProportion', 'var') || isempty(RejectProportion)
  elseif RejectProportion < 0 || RejectProportion > 1
    error('RejectProportion (%%) must be between 0 and 100.')
  else
    set(handles.edit_Proportion, 'String', num2str(RejectProportion, handles.NumFormat));
    handles.Proportion = RejectProportion;
  end
  if ~exist('RejectFitThreshold', 'var') || isempty(RejectFitThreshold)
  elseif RejectFitThreshold < 0 || RejectFitThreshold > 100
    error('RejectFitThreshold (%%) must be between 0 and 100.')
  else
    set(handles.edit_FitThresh, 'String', num2str(RejectFitThreshold, handles.NumFormat));
  end
  
  if ~exist('Renderer', 'var') || isempty(Renderer) 
  elseif ~ismember(Renderer, {'opengl', 'zbuffer', 'painters'})
    error('Renderer must be one of ''opengl'', ''zbuffer'', ''painters'', or empty.');
  elseif strcmpi(Renderer, 'opengl')
    set(hObject, 'Renderer', 'opengl'); % 'opengl', 'zbuffer', 'painters'
    set(handles.Menu_Opengl, 'Checked', 'on');
  else
    set(hObject, 'Renderer', Renderer);
    set(handles.Menu_Opengl, 'Checked', 'off');
  end
  set(hObject, 'DefaultLineLineSmoothing', 'on');
  set(hObject, 'DefaultPatchLineSmoothing', 'on');

  % Apply font size.
  TextObjects = findobj(hObject, '-property', 'FontSize');
  for h = TextObjects
    set(h, 'FontSize', handles.FontSize);
  end
  
  % Specific to Fieldtrip.
  if handles.Fieldtrip
    if ~isfield(handles.ft.cfg, 'dataset') || ~isfield(handles.ft.cfg, 'trl')
      error('Unrecognized Fieldtrip configuration structure.  Should be the output of ft_definetrial');
    end
    Dataset = handles.ft.cfg.dataset;
    
    % Restoring dataset doesn't make sense here, we're just outputting
    % variables.
    set(handles.Menu_Restore, 'Enable', 'off', 'Visible', 'off');
    
    % Extra variable to keep track of rejected trials
    handles.ft.Rejected = [];
  end
  
  % Get data.
  if ~exist('Dataset', 'var') || isempty(Dataset)
    handles.Dataset = '';
    %     handles.Dataset = Menu_Open_Callback(handles.Menu_Open, 0, handles);
  elseif ~ischar(Dataset)
    error('HeadMotionTool:BadDataset', 'Dataset must be a character string, e.g. ''path/dataset.ds''.');
  elseif ~isdir(Dataset)
    error('Dataset not found: %s\n', Dataset);
    %     handles.Dataset = '';
  else
    handles.Dataset = Dataset;
    handles = ReadDataset(handles);
    % If a dataset is provided from the function call, prevent loading
    % another dataset.  This could lead to confusion and errors.  (Wanted
    % to only do this when also output arguments, but can't verify in
    % opening function without a global variable.)
    %     if nargout
    set(handles.Menu_Open, 'Enable', 'off');
    %     end
  end
  
  % ___________________________________

  % Calculation and display will be called either from the
  % button_Correct_Callback or later when a dataset is opened.
  
  % Update handles structure. (Not necessary here but more importantly,
  % cannot be after radio or button calls since they also will update
  % handles.)
  %   guidata(hObject, handles);
  
  if CorrectInitial
    % Check that it's possible, or give error.
    if get(handles.uipanel_Ref, 'SelectedObject') == handles.radio_Initial
      error('Initial position cannot be corrected because it is being used (UseInitial = true).');
    elseif ~handles.CorrectAvailable
      error('Initial position correction is done with CTF command line program changeHeadPos, which seems to have failed (in particular it is only available for UNIX systems).');
    elseif ~isdir(handles.Dataset)
      error('HeadMotionTool:NoDataset', 'Dataset not found for initial position correction.');
    elseif str2double(get(handles.edit_RefDistance, 'String')) < handles.MinDist 
      fprintf('No initial position correction required.\n');
      % Disable correct button then.
      set(handles.button_Correct, 'Enable', 'off');
      % Here we have a valid dataset.
      CalculateMovement(handles);
    else
      % "Press button".
      button_Correct_Callback(hObject, 0, handles);
      % This will also update RefDistance and disable button, and calls
      % CalculateMovement.
    end
  elseif isdir(handles.Dataset)
    % Here we have a valid dataset.
    CalculateMovement(handles);
  end
  % guidata gets updated at end of CalculateMovement, reload.
  handles = guidata(hObject);
  
  % Reject trials but only if saving figure, and only if needed.
  % CalculateMovement above will only have enabled the button if something
  % to reject (thus dataset loaded).  Otherwise button disabled in figure.
  if RejectTrials && strcmpi(get(handles.button_Reject, 'Enable'), 'on')
    if isempty(handles.SavePictureFile)
      error('Cannot automatically reject trials without saving figure screenshot for later quality check.  Please provide SavePictureFile.');
    end
    % "Press button".
    button_Reject_Callback(hObject, 0, handles);
  end
  
  % Used to close figure here if ~handles.ShowGUI, but this would mean no
  % output...  Not sure if anybody ever used this seemingly buggy feature
  % (though possibly a change in Matlab version?).
  if handles.ShowGUI 
    % This makes the current execution thread (not the whole gui) halt and
    % wait for the uiresume command.  After the end of this function, the
    % ouput function eventually runs, and execution returns to the calling
    % thread.  Could make the figure modal since command window seems
    % available but is really stuck in "busy".  But unnecessary and can make
    % it harder to deal with potential figure closing function bugs.
    %   hObject.WindowStyle = 'modal';
    uiwait(hObject);
    % However, this means we're "stuck" in the initialization phase of
    % gui_mainfcn, so there's room for improvement.
  % else the GUI will show up on screen to allow a screenshot, but it will
  % close right after, from the output function.  If we don't output nor
  % save a screenshot, we could close it here, but can't easily know if we
  % have outputs.
  end
    
  
end


% --- Executes when user attempts to close figure_HeadMotion.
function figure_HeadMotion_CloseRequestFcn(hObject, eventdata, handles)
  % hObject    handle to figure_HeadMotion (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  % Resume execution of main thread but in parallel or after finishing this
  % function.  (Thus if we deleted the figure here after resuming, it would
  % happen before the output function, which would then have an emtpy
  % handles structure.)  
  if isequal(get(hObject, 'waitstatus'), 'waiting')
    %  The GUI is "in UIWAIT", so "UIRESUME"
    % As uiwait is in the opening function, the output function will get
    % called after this, and then return execution to the calling program.
    uiresume(hObject);
  else
    % Hint: delete(hObject) closes the figure
    delete(hObject);
  end
  
  % No longer save user parameters here because of potential issues running
  % this on a server with multiple users.  Must now be saved manually with
  % menu.
    
end


% --- Outputs from this function are returned to the command line.
function varargout = HeadMotionTool_OutputFcn(hObject, eventdata, handles)
  % varargout  cell array for returning output args (see VARARGOUT);
  % hObject    handle to figure
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  % Save screenshot of GUI if wanted.
  if ~isempty(handles.SavePictureFile)
    if exist(handles.SavePictureFile, 'file')
      warning('Image file already exists, adding timestamp to file name.');
      [Path, Name] = fileparts(handles.SavePictureFile);
      handles.SavePictureFile = fullfile(Path, [Name, '_', datestr(now, 30), '.png']);
    end
    % Matlab is horrible at saving figures, shame on it!
    % Problems saving the figure:
    % 1. saveas: Resizes.
    % 2. print: Missing text if figure is not visible.
    % 3. print: Bad looking text (no font smoothing) if opengl.  painters
    %    ok but thinner fonts than on screen.
    % 4. export_fig: either aa everything, making crosshair look faint
    %    and transparent, or same as print (opengl or painters issues
    %    depending on whether there is an overlay).
    % 5. screencapture: need figure border.
    
    %       set(hObject, 'Visible', 'on'); % Already visible here.
    Pos = get(hObject, 'Position');
    if numel(handles.WindowFrameSize) == 1
      handles.WindowFrameSize = handles.WindowFrameSize * [1, 1];
    end
%     screencapture(hObject, [handles.WindowFrameSize, Pos([3, 4])], ...
%       handles.SavePictureFile);
    saveas(hObject,handles.SavePictureFile);
  end

  % Output if a dataset was loaded.
  if isfield(handles, 'MedianLoc')
    switch get(handles.uipanel_Ref, 'SelectedObject')
      case handles.radio_Median
        varargout{1} = handles.MedianLoc;
      case handles.radio_Initial
        varargout{1} = handles.Fiducials;
    end
    varargout{2} = str2double(get(handles.edit_MaxDistance, 'String'));
    % Figure is now modal, it is gone once we return outputs, no need for
    % figure handle.
    %     varargout{3} = hObject;
    if handles.Fieldtrip
      % Remove rejected trials.
      handles.ft.cfg.trl(handles.ft.Rejected, :) = [];
      varargout{3} = handles.ft.cfg;
      varargout{4} = handles.ft.grad;
    elseif nargout > 2
      error('cfg and grad outputs are only available with Fieldtrip data, provide the ''Fieldtrip'' input parameter on the command line or remove the extra output variables when calling the function.');
    end
  elseif nargout > 0
    error('The program can only return output if all required arguments (at least a valid dataset or Fieldtrip "definetrial" structure) are provided as inputs on the command line.');
  end

  % Otherwise, don't return anything.

  % We close the figure here because after resuming from "uiwait", if it
  % deletes right away, handles is gone once we get to the output function.
  close(hObject);

end


% --- Executes when figure_HeadMotion is resized.
function figure_HeadMotion_ResizeFcn(hObject, eventdata, handles)
  % hObject    handle to figure_HeadMotion (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  
  % Can get called before it should, so add safety.
  if ~isfield(handles, 'PanelPositions')
    return;
  end
  
  % ---------------------------------------------------------------------
  % Play with positions.
  
  % Get new figure size.
  FigPos = get(hObject, 'Position');
  
  % Objects at top only move vertically.
  for h = [handles.togglebutton_InterCoil, handles.uipanel_InterCoil, ...
      handles.uipanel_Thresh, handles.uipanel_Global]
    Pos = get(h, 'Position');
    Pos(2) = FigPos(4) - handles.PanelPositions(1);
    set(h, 'Position', Pos);
  end
  h = handles.uipanel_Ref;
  Pos = get(h, 'Position');
  Pos(2) = FigPos(4) - handles.PanelPositions(2);
  set(h, 'Position', Pos);
  
  % Axes need to be stretched.
  for h = [handles.axes_Trials, handles.axes_TrialsFit, ...
      handles.axes_InterCoil]
    Pos = get(h, 'Position');
    Pos(3) = FigPos(3) - 2 * handles.DefaultAxisPosition(1);
    Pos(4) = FigPos(4) - handles.DefaultAxisPosition(2) - handles.AxisTop;
    set(h, 'Position', Pos);
  end
  
  % Axes labels.
  if isfield(handles, 'YLabel')
    for h = [handles.YLabel, handles.YLabelInterCoil, handles.YLabelFit]
      TxtPos = get(h, 'Position');
      TxtPos(2) = Pos(4)/2;
      set(h, 'Position', TxtPos);
    end
    % YLabelFit also needs to move in x.
    TxtPos(1) = Pos(3) + handles.TextPos;
    set(h, 'Position', TxtPos);
    for h = [handles.XLabel, handles.XLabelInterCoil]
      TxtPos = get(h, 'Position');
      TxtPos(1) = Pos(3)/2;
      set(h, 'Position', TxtPos);
    end
  end
  
  % Slider width.
  Pos = get(handles.slider_Trials, 'Position');
  Pos(3) = FigPos(3) - Pos(1) - handles.DefaultAxisPosition(1);
  set(handles.slider_Trials, 'Position', Pos);
end





% -----------------------------------------------------------------------
% Empty object functions.
% -----------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function slider_Trials_CreateFcn(hObject, eventdata, handles)
  % hObject    handle to slider_Trials (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    empty - handles not created until after all CreateFcns called
  
  % Hint: slider controls usually have a light gray background.
  if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
  end
end

% --- Executes during object creation, after setting all properties.
function edit_MaxFit_CreateFcn(hObject, eventdata, handles)
  % hObject    handle to edit_MaxFit (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    empty - handles not created until after all CreateFcns called
  
  % Hint: edit controls usually have a white background on Windows.
  %       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
end

% --- Executes during object creation, after setting all properties.
function edit_Distance_Q_CreateFcn(hObject, eventdata, handles)
  % hObject    handle to edit_Distance_Q (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    empty - handles not created until after all CreateFcns called
  
  % Hint: edit controls usually have a white background on Windows.
  %       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
end

% --- Executes during object creation, after setting all properties.
function edit_MaxDistance_CreateFcn(hObject, eventdata, handles)
  % hObject    handle to edit_MaxDistance (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    empty - handles not created until after all CreateFcns called
  
  % Hint: edit controls usually have a white background on Windows.
  %       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
end

% --- Executes during object creation, after setting all properties.
function edit_Range_CreateFcn(hObject, eventdata, handles)
  % hObject    handle to edit_Range (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    empty - handles not created until after all CreateFcns called
  
  % Hint: edit controls usually have a white background on Windows.
  %       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
end

% --- Executes during object creation, after setting all properties.
function edit_FitThresh_CreateFcn(hObject, eventdata, handles)
  % hObject    handle to edit_FitThresh (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    empty - handles not created until after all CreateFcns called
  
  % Hint: edit controls usually have a white background on Windows.
  %       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
end

% --- Executes during object creation, after setting all properties.
function edit_Proportion_CreateFcn(hObject, eventdata, handles)
  % hObject    handle to edit_Proportion (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    empty - handles not created until after all CreateFcns called
  
  % Hint: edit controls usually have a white background on Windows.
  %       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
end

% --- Executes during object creation, after setting all properties.
function edit_Thresh_CreateFcn(hObject, eventdata, handles)
  % hObject    handle to edit_Thresh (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    empty - handles not created until after all CreateFcns called
  
  % Hint: edit controls usually have a white background on Windows.
  %       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
end

% --- Executes during object creation, after setting all properties.
function edit_RefDistance_CreateFcn(hObject, eventdata, handles)
  % hObject    handle to edit_RefDistance (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    empty - handles not created until after all CreateFcns called
  
  % Hint: edit controls usually have a white background on Windows.
  %       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
end

% --- This function must be there even if empty.
function Menu_Callback(hObject, eventdata, handles)
  % hObject    handle to Menu (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
end


% --- Executes during object creation, after setting all properties.
function edit_nI_CreateFcn(hObject, eventdata, handles)
  % hObject    handle to edit_nI (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    empty - handles not created until after all CreateFcns called
  
  % Hint: edit controls usually have a white background on Windows.
  %       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
end


% --- Executes during object creation, after setting all properties.
function edit_InterCoilThresh_CreateFcn(hObject, eventdata, handles)
  % hObject    handle to edit_InterCoilThresh (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    empty - handles not created until after all CreateFcns called
  
  % Hint: edit controls usually have a white background on Windows.
  %       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
end


% --- Executes during object creation, after setting all properties.
function edit_InterCoilTrial_CreateFcn(hObject, eventdata, handles)
  % hObject    handle to edit_InterCoilTrial (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    empty - handles not created until after all CreateFcns called
  
  % Hint: edit controls usually have a white background on Windows.
  %       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
end


% --- Executes during object creation, after setting all properties.
function edit_InterCoilSample_CreateFcn(hObject, eventdata, handles)
  % hObject    handle to edit_InterCoilSample (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    empty - handles not created until after all CreateFcns called
  
  % Hint: edit controls usually have a white background on Windows.
  %       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
end


function edit_NALE_Callback(hObject, eventdata, handles)
  % hObject    handle to edit_NALE (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  % Hints: get(hObject,'String') returns contents of edit_NALE as text
  %        str2double(get(hObject,'String')) returns contents of edit_NALE as a double
end


% --- Executes during object creation, after setting all properties.
function edit_NALE_CreateFcn(hObject, eventdata, handles)
  % hObject    handle to edit_NALE (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    empty - handles not created until after all CreateFcns called
  
  % Hint: edit controls usually have a white background on Windows.
  %       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
end


function edit_LERE_Callback(hObject, eventdata, handles)
  % hObject    handle to edit_LERE (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  % Hints: get(hObject,'String') returns contents of edit_LERE as text
  %        str2double(get(hObject,'String')) returns contents of edit_LERE as a double
end


% --- Executes during object creation, after setting all properties.
function edit_LERE_CreateFcn(hObject, eventdata, handles)
  % hObject    handle to edit_LERE (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    empty - handles not created until after all CreateFcns called
  
  % Hint: edit controls usually have a white background on Windows.
  %       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
end


function edit_RENA_Callback(hObject, eventdata, handles)
  % hObject    handle to edit_RENA (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  % Hints: get(hObject,'String') returns contents of edit_RENA as text
  %        str2double(get(hObject,'String')) returns contents of edit_RENA as a double
end

% --- Executes during object creation, after setting all properties.
function edit_RENA_CreateFcn(hObject, eventdata, handles)
  % hObject    handle to edit_RENA (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    empty - handles not created until after all CreateFcns called
  
  % Hint: edit controls usually have a white background on Windows.
  %       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
end




% -----------------------------------------------------------------------
% Object callback functions.
% -----------------------------------------------------------------------


% --- Executes when selected object is changed in uipanel_Ref.
function uipanel_Ref_SelectionChangeFcn(hObject, eventdata, handles)
  % hObject    handle to the selected object in uipanel_Ref
  % eventdata  structure with the following fields (see UIBUTTONGROUP)
  %	EventName: string 'SelectionChanged' (read only)
  %	OldValue: handle of the previously selected object or empty if none was selected
  %	NewValue: handle of the currently selected object
  % handles    structure with handles and user data (see GUIDATA)
  
  % RefDistance is initially 0, so this will ensure we do nothing until a
  % dataset is loaded.  Also, if a dataset is loaded but the value is 0,
  % nothing changes when we choose Median or Initial.
  if str2double(get(handles.edit_RefDistance, 'String')) < handles.MinDist  
    % Nothing to do.  (Exclusive selection is automatic.)
    return
  end
  
  switch get(eventdata.NewValue, 'Tag') % Get Tag of selected object.
    case 'radio_Median'
      if handles.CorrectAvailable
        % Enable button.
        set(handles.button_Correct, 'Enable', 'on');
      end
    case 'radio_Initial'
      % Disable button.
      set(handles.button_Correct, 'Enable', 'off');
    otherwise
      error('HeadMotionTool:BadRadio', 'Unrecognized selection for reference position.');
  end
  
  % Need to explicitely disable radio buttons during the following
  % calculations, or they will behave like check boxes (independently
  % (un)selectable), leading to confusion and errors.
  set(handles.radio_Initial, 'Enable', 'off');
  set(handles.radio_Median, 'Enable', 'off');
  CalculateMovement(handles);
  set(handles.radio_Initial, 'Enable', 'on');
  set(handles.radio_Median, 'Enable', 'on');
end
  
  
% --- Executes on button press in button_Correct.
function button_Correct_Callback(hObject, eventdata, handles)
  % hObject    handle to button_Correct (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  if handles.Fieldtrip
    % Correct sensor positions, which are in "initial" CTF head coordinates
    % in the grad structure.
    % Reference fiducial locations in "initial" CTF head coordinates, cm.
    NewFids = ChangeCoordinates(reshape(handles.MedianLoc, [3, 3])', ...
      reshape(handles.Fiducials, [3, 3])', 0.1);
    
    for f = {'chanpos', 'coilpos'}
      handles.ft.grad.(f{1}) = ChangeCoordinates(handles.ft.grad.(f{1}), NewFids);
    end
    for f = {'chanori', 'coilori'}
      handles.ft.grad.(f{1}) = ChangeCoordinates(handles.ft.grad.(f{1}), NewFids, ...
        1, 1, 'CTF', false, true, false);
        %     Scaling, Orientation, System, Inverse, Vector, Verbose)
    end
    
    handles.Fiducials = handles.MedianLoc;
        
  else
    % Make sure there is a backup of the original .hc file.
    [unused, DatasetName] = fileparts(handles.Dataset);
    BackupOriginal([handles.Dataset, filesep, DatasetName, '.hc']);
    
    % Correct initial position with CTF command line program changeHeadPos.
    % changeHeadPos -na x y z -le x y z -re x y z dataset.ds  (dewar coordinates in mm)
    CommandString = sprintf('changeHeadPos -na %f %f %f -le %f %f %f -re %f %f %f %s', ...
      handles.MedianLoc, handles.Dataset);
    [Status, Output] = system(CommandString);
    if Status
      fprintf(Output);
      error('Initial position correction error.');
    else
      fprintf('Initial position correction completed successfully.\n');
    end
    
    % Update initial position to be the same as the reference, otherwise if
    % we reject trials after, the old initial position would be reused.
    handles.Fiducials = openhc(handles.Dataset, true); % 3x3 in cm, dewar coordinates.
    % Convert to mm and one line.
    handles.Fiducials = 10 * reshape(handles.Fiducials', 1, []);
  end
  
  % Update RefDistance and disable button.
  set(handles.edit_RefDistance, 'String', '0.00');
  set(handles.button_Correct, 'Enable', 'off');

  CalculateMovement(handles); % updates guidata at end
  
end


function edit_FitThresh_Callback(hObject, eventdata, handles)
  % hObject    handle to edit_FitThresh (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  NewValue = str2double(get(hObject, 'String'));
  if isnan(NewValue) || isempty(NewValue) || NewValue < 0 || NewValue > 100
    beep;
    fprintf('Fit error (%%) threshold must be between 0 and 100.\n');
    % Reset to old value, from threshold line if it exists or default.
    if isfield(handles, 'Plot_FitThresh')
      NewValue = get(handles.Plot_FitThresh, 'YData');
      set(hObject, 'String', num2str(NewValue(1), handles.NumFormat));
    else
      set(hObject, 'String', '10.0');
    end
    return;
  else
    % Apply our format.
    set(hObject, 'String', num2str(NewValue, handles.NumFormat));
  end

  % Update everything if a dataset is loaded.
  if isfield(handles, 'Plot_LocationDist')
    CalculateMovement(handles);
  end

end


function edit_Thresh_Callback(hObject, eventdata, handles)
  % hObject    handle to edit_Thresh (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  NewValue = str2double(get(hObject, 'String'));
  if isnan(NewValue) || isempty(NewValue) || NewValue < 0
    beep;
    fprintf('Distance (mm) threshold must be a positive number.\n');
    % Reset to old value, from threshold line if it exists or default.
    if isfield(handles, 'Plot_Thresh')
      NewValue = get(handles.Plot_Thresh, 'YData');
      set(hObject, 'String', num2str(NewValue(1), handles.NumFormat));
    else
      set(hObject, 'String', '5.0');
    end
    return;
  else
    % Apply our format.
    set(hObject, 'String', num2str(NewValue, handles.NumFormat));
  end

  % Update everything if a dataset is loaded.
  if isfield(handles, 'Plot_LocationDist')
    CalculateMovement(handles);
  end

end


function edit_Proportion_Callback(hObject, eventdata, handles)
  % hObject    handle to edit_Proportion (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  NewValue = str2double(get(hObject, 'String'));
  if isnan(NewValue) || isempty(NewValue) || NewValue < 0 || NewValue > 100
    beep;
    fprintf('Proportion (%%) must be between 0 and 100.\n');
    % Reset to old value.
    set(hObject, 'String', num2str(handles.Proportion, handles.NumFormat));
    return;
  else
    % Apply our format.
    set(hObject, 'String', num2str(NewValue, handles.NumFormat));
  end
  handles.Proportion = NewValue;

  % Update everything if a dataset is loaded.
  if isfield(handles, 'Plot_LocationDist')
    CalculateMovement(handles);
  end

end


% --- Executes on button press in button_Reject.
function button_Reject_Callback(hObject, eventdata, handles)
  % hObject    handle to button_Reject (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  % ~isempty(NewBadTrials) required for button to be active.
  
  if handles.Fieldtrip
    % Remove unwanted trials from Fieldtrip structure. Only keep track of
    % them for now and remove when returning output, otherwise keeping
    % track of trial numbers with repeated rejection steps would be
    % difficult.
    for bad = 1:numel(handles.NewClass)
      handles.ft.Rejected = [handles.ft.Rejected; handles.NewClass(bad).Trials];
    end
  else
    % Save bad trial information in ClassFile.cls in dataset.
    % Could use CTF command line program: "headMotionDetect -fe 30 -pc 10
    % -uf dataset.ds", but it's only a matter of writing "bad" trial
    % classes in ClassFile.cls.
    savecls(handles.NewClass, handles.Dataset);
  end
  
  handles.Class = handles.NewClass;
  
  % Update list of good trials, which is used to find new vs. old bad
  % trials.
  handles.GoodTrials = GetGoodTrials(handles.Class, handles.nT);
  
  % [THIS IS NO LONGER THE CASE]
  % Get rigid reference body, in head coordinates.  We do it again here to
  % get consistent trials/samples in I, and also the Rigid body may change
  % if we rejected initial trials.
  %handles = ReferenceBody(handles);
  
  % Get median location for each coil separately.
  handles.MedianLoc = MedianLocation(handles.Data, handles.GoodTrials);
  % Correct the median coil locations to impose rigid body shape.
  handles.MedianLoc(:) = CorrectRigid(handles.MedianLoc, handles.Rigid);
  % Old way, with full ChangeCoordinates function:
  %   handles.MedianLoc(:) = ChangeCoordinates(...
  %     reshape(handles.Rigid, [3, 3])', reshape(handles.MedianLoc, [3, 3])', ...
  %     [], [], [], true)'; % Inverse is true.

  % Update reference position objects.
  if strcmpi(get(handles.Menu_UseRigid, 'Checked'), 'on') % handles.UseRigid
    D = RigidDistances(handles.MedianLoc, handles.Fiducials);
  else
    D = MaxCoilNorms(handles.MedianLoc - handles.Fiducials);
  end
  set(handles.edit_RefDistance, 'String', num2str(D , handles.NumFormatSmall));
  if D > handles.MinDist && ... 
      get(handles.uipanel_Ref, 'SelectedObject') == handles.radio_Median && ...
      handles.CorrectAvailable
    set(handles.button_Correct, 'Enable', 'on');
  else
    set(handles.button_Correct, 'Enable', 'off');
  end
  
  CalculateMovement(handles); % updates guidata at end
  
end


% --------------------------------------------------------------------
function Menu_Open_Callback(hObject, eventdata, handles)
  % hObject    handle to Menu_Open (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  %   fprintf([handles.Dataset, '\n']);
  % Get dataset
  Dataset = uigetdir(handles.DatasetDir, 'Select dataset.');
  if ~ischar(Dataset) || ~isdir(Dataset)
    % User probably hit cancel, do nothing.
    return;
  end
  
  handles.Dataset = Dataset;
  handles = ReadDataset(handles);
    
  % Calculate and update everything else (Graphs, overall values and
  % state of Reject button).
  %   if length(handles.GoodTrials) > 1
  CalculateMovement(handles);
  %   end

end


% --------------------------------------------------------------------
function Menu_Restore_Callback(hObject, eventdata, handles)
  % hObject    handle to Menu_Restore (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  Confirmation = questdlg('Restore original dataset initial position and trial classification?', ...
    'Restore dataset...', 'Ok', 'Cancel', 'Ok');
  
  if strcmpi(Confirmation, 'Ok')
    % opencls(ds, true) restores the original ClassFile file if present.
    opencls(handles.Dataset, true, false); % Quiet
    
    [Path, Name] = fileparts(handles.Dataset);
    HCFileOriginal = [handles.Dataset, filesep, Name, '_Original.hc'];
    if handles.CorrectAvailable && exist(HCFileOriginal, 'file');
      % openhc(..., true) does not restore the file since sensor coordinates
      % in the res4 file need to be modified as well.
      Fiducials = openhc(handles.Dataset, true, [], true); % 3x3 in cm, dewar coordinates, original.
      % Convert to mm and one line.
      Fiducials = 10 * reshape(Fiducials', 1, []);
      % Correct initial position with CTF command line program changeHeadPos.
      % changeHeadPos -na x y z -le x y z -re x y z dataset.ds  (dewar coordinates in mm)
      CommandString = sprintf('changeHeadPos -na %f %f %f -le %f %f %f -re %f %f %f %s', ...
        Fiducials, handles.Dataset);
      [Status, Output] = system(CommandString);
      if Status
        fprintf(Output);
        error('Initial position restoration error.');
      else
        fprintf('Initial position restoration completed successfully.\n');
      end
    end
    
    % Reload dataset.
    handles = ReadDataset(handles);
    CalculateMovement(handles);
  end
      
end


% --------------------------------------------------------------------
function Menu_Toolbar_Callback(hObject, eventdata, handles)
  % hObject    handle to Menu_Toolbar (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  if strcmpi(get(hObject, 'Checked'), 'on')
    % Turn off.
    set(hObject, 'Checked', 'off');
    set(handles.figure_HeadMotion, 'MenuBar', 'none');
    set(handles.figure_HeadMotion, 'ToolBar', 'auto');
  else
    % Turn on.
    set(hObject, 'Checked', 'on');
    set(handles.figure_HeadMotion, 'MenuBar', 'figure');
    set(handles.figure_HeadMotion, 'ToolBar', 'figure');
  end
end


% --------------------------------------------------------------------
function Menu_Opengl_Callback(hObject, eventdata, handles)
  % hObject    handle to Menu_Opengl (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  if strcmpi(get(hObject, 'Checked'), 'on')
    % Turn off.
    set(hObject, 'Checked', 'off');
    set(handles.figure_HeadMotion, 'Renderer', 'painters');
  else
    % Turn on.
    set(hObject, 'Checked', 'on');
    set(handles.figure_HeadMotion, 'Renderer', 'opengl');
  end
end


% --------------------------------------------------------------------
function Menu_FontSize_Callback(hObject, eventdata, handles)
  % hObject    handle to Menu_FontSize (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  FontSize = inputdlg(['Desired font size (currently ', ...
    num2str(handles.FontSize, '%0.1f'), ' pts):']);
  if isempty(FontSize)
    % Cancelled.
    return;
  end
  
  FontSize = str2double(FontSize{1});
  if isnan(FontSize) || ~isnumeric(FontSize) || FontSize < 1
    % Bad user input; do nothing.
    return;
  end
  
  handles.FontSize = FontSize;
  
  for h = (findobj(handles.figure_HeadMotion, '-property', 'FontSize'))'
    set(h, 'FontSize', handles.FontSize);
  end

  handles.TextOptions{2} = handles.FontSize;

  % Save updated handles data.
  guidata(hObject, handles);
end


% --------------------------------------------------------------------
function Menu_UseRigid_Callback(hObject, eventdata, handles)
  % hObject    handle to Menu_UseRigid (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  if strcmpi(get(hObject, 'Checked'), 'on')
    % Turn off.
    set(hObject, 'Checked', 'off');
  else
    % Turn on.
    set(hObject, 'Checked', 'on');
  end
  % For simplicity, restart from beginning: read dataset and update.
  handles = ReadDataset(handles);
  CalculateMovement(handles);
end


% --------------------------------------------------------------------
function Menu_Help_Callback(hObject, eventdata, handles)
  % hObject    handle to Menu_Help (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  helpdlg(['For (limited) help on this tool, type "help HeadMotionTool" on the Matlab command line.  ', ...
    'For additional help, please contact Marc Lalancette (x1535).  ', ...
    '(A more detailed help document could be created if there is demand for it.)'], ...
    'HeadMotionTool Help');
end


% --------------------------------------------------------------------
function Menu_Exit_Callback(hObject, eventdata, handles)
  % hObject    handle to Menu_Exit (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  % This just calls the figure's CloseRequestFcn.
  close(handles.figure_HeadMotion);
end


% --- Executes on slider movement.
function slider_Trials_Callback(hObject, eventdata, handles) %#ok<*DEFNU,*INUSD>
  % hObject    handle to slider_Trials (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  Value = round(get(hObject, 'Value'));
  % Adjust X limits  on both axes based on slider range and value.
  XLim = Value + [0, ...
    handles.nT - get(hObject, 'Max') + 1];
  set(hObject, 'Value', Value);
  for h = [handles.axes_Trials, handles.axes_TrialsFit, ...
      handles.axes_InterCoil]
    set(h, 'XLim', XLim);
  end
end


% --- Executes on button press in pushbutton_ZoomIn.
function pushbutton_ZoomIn_Callback(hObject, eventdata, handles)
  % hObject    handle to pushbutton_ZoomIn (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  Max = get(handles.slider_Trials, 'Max');
  
  XWidth = round((1 + handles.nT - Max) / handles.Phi); % This gives a
  % minimum of 1 when applied repeatedly to any number >= 1.
  Max = 1 + handles.nT - XWidth;
  % When zooming in, Value doesn't need to change.
  %   Value = min(Max, get(handles.slider_Trials, 'Value'));
  
  % Adjust slider range and value.
  set(handles.slider_Trials, 'Max', Max);
  if Max ~= 1
    set(handles.slider_Trials, 'Enable', 'on');
    set(handles.slider_Trials, 'Min', 1);
  end
  set( handles.slider_Trials, 'SliderStep', ...
    [1, XWidth]/(Max - get(handles.slider_Trials, 'Min')) ); 
  %   set(handles.slider_Trials, 'Value', Value);
  
  slider_Trials_Callback(handles.slider_Trials, 0, handles)
  
end


% --- Executes on button press in pushbutton_ZoomOut.
function pushbutton_ZoomOut_Callback(hObject, eventdata, handles)
  % hObject    handle to pushbutton_ZoomOut (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  Max = get(handles.slider_Trials, 'Max');
  
  XWidth = min(handles.nT, round((1 + handles.nT - Max) * handles.Phi)); % This gives the
  % Fibonacci numbers when applied repeatedly starting with 1.
  Max = 1 + handles.nT - XWidth;
  Value = min(Max, get(handles.slider_Trials, 'Value'));
  
  % Adjust slider range and value.
  if Max == 1
    set(handles.slider_Trials, 'Enable', 'off');
    set(handles.slider_Trials, 'Min', 0);
  end
  set(handles.slider_Trials, 'Max', Max);
  set( handles.slider_Trials, 'SliderStep', ...
    [1, XWidth]/(Max - get(handles.slider_Trials, 'Min')) ); 
  set(handles.slider_Trials, 'Value', Value);

  slider_Trials_Callback(handles.slider_Trials, 0, handles)
  
end


% --- Executes on button press in togglebutton_InterCoil.
function togglebutton_InterCoil_Callback(hObject, eventdata, handles)
  % hObject    handle to togglebutton_InterCoil (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  % Hint: get(hObject,'Value') returns toggle state of togglebutton_InterCoil
  
  if get(hObject, 'Value')
    % Button was pressed.
    % Show axis and pannel.
    for h = [handles.axes_InterCoil, handles.uipanel_InterCoil] %, handles.XLabelInterCoil, handles.YLabelInterCoil
      set(h , 'Visible', 'on');
      set(get(h, 'children'), 'Visible', 'on');
    end
    % Hide other axes and panels.
    for h = [handles.axes_Trials, handles.axes_TrialsFit, ...
        handles.uipanel_Thresh, handles.uipanel_Global]  %, handles.XLabel, handles.YLabel, handles.YLabelFit
      set(h, 'Visible', 'off');
      set(get(h, 'children'), 'Visible', 'off');
    end
    set(hObject, 'String', 'Back to head motion');
    % Except bad trial highlighting, keep those visible.
    if isfield(handles, 'BadTrPatches')
      set(handles.BadTrPatches, 'Visible', 'on');
    end
    if isfield(handles, 'NewBadTrPatches')
      set(handles.NewBadTrPatches, 'Visible', 'on');
    end
    % And axes_TrialsFit itself to keep white background, but need to "hide" axis.
    set(handles.axes_TrialsFit, 'Visible', 'on', 'YTick', []);
    %     axis(handles.axes_TrialsFit, 'off'); % This seems to be the same
                                % property as above, i.e. axes, not axis.
  else
    % Button was released.
    % Hide axis and pannel.
    for h = [handles.axes_InterCoil, handles.uipanel_InterCoil] %, handles.XLabelInterCoil, handles.YLabelInterCoil
      set(h , 'Visible', 'off');
      set(get(h, 'children'), 'Visible', 'off');
    end
    % Show other axes and panels.
    for h = [handles.axes_Trials, handles.axes_TrialsFit, ...
        handles.uipanel_Thresh, handles.uipanel_Global] %, handles.XLabel, handles.YLabel, handles.YLabelFit
      set(h, 'Visible', 'on');
      set(get(h, 'children'), 'Visible', 'on');
    end
    set(hObject, 'String', 'Verify coils');
    % Need to reset y axis on TrialsFit.
    set(handles.axes_TrialsFit, 'YTickMode', 'auto');
  end
  
end


function edit_InterCoilThresh_Callback(hObject, eventdata, handles)
  % hObject    handle to edit_Thresh (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  NewValue = str2double(get(hObject, 'String'));
  if isnan(NewValue) || isempty(NewValue) || NewValue < 0
    beep;
    fprintf('Distance (mm) threshold must be a positive number.\n');
    % Reset to default.
    set(hObject, 'String', '0.5');
    return;
  else
    % Apply our format.
    set(hObject, 'String', num2str(NewValue, handles.NumFormat));
  end

  % Update everything if a dataset is loaded.
  if isfield(handles, 'Plot_LocationDist')
    handles = ReferenceBody(handles);
    CalculateMovement(handles);
  end

end


% --------------------------------------------------------------------
function Menu_LoadParam_Callback(hObject, eventdata, handles)
  % hObject    handle to Menu_LoadParam (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)

  [File, Path] = uigetfile(handles.UserDataFile, 'Load parameters');
  if File == 0
    return;
  end
  handles.UserDataFile = fullfile(Path, File);
 
  load(handles.UserDataFile);
  % Load user parameters in text boxes.
  for EditBox = (fieldnames(UserData.Edit))'
    set(handles.(EditBox{1}), 'String', UserData.Edit.(EditBox{1}));
  end
  
  if UserData.uipanel_Ref
    set(handles.uipanel_Ref, 'SelectedObject', handles.radio_Initial);
  else
    set(handles.uipanel_Ref, 'SelectedObject', handles.radio_Median);
  end
  
  % Get directory of last loaded dataset.
  handles.DatasetDir = UserData.DatasetDir;
  % Other options.
  handles.FontSize = UserData.FontSize;
  for h = (findobj(handles.figure_HeadMotion, '-property', 'FontSize'))'
    set(h, 'FontSize', handles.FontSize);
  end
  handles.TextOptions{2} = handles.FontSize;
  
  set(handles.Menu_Opengl, 'Checked', UserData.Menu_Opengl);
  if strcmpi(UserData.Menu_Opengl, 'on')
    set(handles.figure_HeadMotion, 'Renderer', 'opengl');
  else
    set(handles.figure_HeadMotion, 'Renderer', 'painters');
  end
  
  % Save updated handles data.
  guidata(hObject, handles);
end


% --------------------------------------------------------------------
function Menu_SaveParam_Callback(hObject, eventdata, handles)
  % hObject    handle to Menu_SaveParam (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)

  [File, Path] = uigetfile(handles.UserDataFile, 'Save parameters');
  if File == 0
    return;
  end
  handles.UserDataFile = fullfile(Path, File);
 
  % Save user parameters.
  %   if handles.ShowGUI
  for h = (findobj(handles.figure_HeadMotion, 'Style', 'edit', 'Enable', 'on'))'
    UserData.Edit.(get(h, 'Tag')) = get(h, 'String');
  end
  % true means Initial is selected, false means Median.
  UserData.uipanel_Ref = get(handles.uipanel_Ref, 'SelectedObject') == handles.radio_Initial;
  % Save directory of last loaded dataset.
  UserData.DatasetDir = handles.DatasetDir; % fileparts(handles.Dataset);
  % Other options.
  UserData.FontSize = handles.FontSize;
  UserData.Menu_Opengl = get(handles.Menu_Opengl, 'Checked');
  
  % Save the mat file in the same directory as the program (this m-file).
  save(handles.UserDataFile, 'UserData');

end



% -----------------------------------------------------------------------
% Other functions.
% -----------------------------------------------------------------------


% ----------------------------------------------------------------------
% Wrapper for GetData function, to also update reference position objects.
function handles = ReadDataset(handles)
  %   handles.Dataset = Dataset; % Used to be second input argument.
  [handles.DatasetDir, DatasetName] = fileparts(handles.Dataset);
  set(handles.figure_HeadMotion, 'Name', ['Head Motion Tool - ', DatasetName]);
  
  %   [handles.Fiducials, ...
  %     handles.Data, handles.FitData, handles.nS, handles.nT, ...
  %     handles.Class, handles.GoodTrials] = ...
  %     GetData(handles.Dataset);
  handles = GetData(handles);
  % Everything in dewar coordinates, in mm.
  % MedianLoc, Fids: [1, nC, nT]
  % Data: [nS, nC, nT]
  
  % Update scroll bar.
  set(handles.slider_Trials, 'Min', 0);
  set(handles.slider_Trials, 'Max', 1);
  set(handles.slider_Trials, 'Enable', 'off');
  set(handles.slider_Trials, 'SliderStep', [1/handles.nT, 1]); 
  set(handles.slider_Trials, 'Value', 1);

  % Get rigid reference body, in head coordinates.
  pause(1);
  handles = ReferenceBody(handles);
  % Get median location for each coil separately.
  handles.MedianLoc = MedianLocation(handles.Data, handles.GoodTrials);
  % Correct the median coil locations to impose rigid body shape.
  handles.MedianLoc(:) = CorrectRigid(handles.MedianLoc, handles.Rigid);
  % Old way, with full ChangeCoordinates function:
  %   handles.MedianLoc(:) = ChangeCoordinates(...
  %     reshape(handles.Rigid, [3, 3])', reshape(handles.MedianLoc, [3, 3])', ...
  %     [], [], [], true)'; % Inverse is true.
  
  % Update reference position objects.
  if strcmpi(get(handles.Menu_UseRigid, 'Checked'), 'on') % handles.UseRigid
    D = RigidDistances(handles.MedianLoc, handles.Fiducials);
  else
    D = MaxCoilNorms(handles.MedianLoc - handles.Fiducials);
  end
  set(handles.edit_RefDistance, 'String', num2str(D , handles.NumFormatSmall));
  if D > handles.MinDist && ... 
      get(handles.uipanel_Ref, 'SelectedObject') == handles.radio_Median && ...
      handles.CorrectAvailable
    set(handles.button_Correct, 'Enable', 'on');
  else
    set(handles.button_Correct, 'Enable', 'off');
  end
  
end


% ----------------------------------------------------------------------
% Get additional information on movement overall and per trial, possibly
% for automated trial rejection, or to get values that can be used in
% correlations or as weights in group averaging for example.
function CalculateMovement(handles)
  % Set mouse pointer as busy.
  set(handles.figure_HeadMotion, 'Pointer', 'watch');
  drawnow;
  
  nT = handles.nT;
  switch get(handles.uipanel_Ref, 'SelectedObject')
    case handles.radio_Median
      Reference = handles.MedianLoc;
    case handles.radio_Initial
      Reference = handles.Fiducials;
  end

  % Calculate distance from chosen location (median or fids) across trials
  % and time.  Use geometric median location per trial for unique trial
  % distance. Note that this is not an actual position visited in the
  % trial, so it could be outside the range of MinDistance to MaxDistance.
  T_Location = GeoMedian(reshape(handles.Data, [handles.nS, 3, 3, nT]), 1e-3); % [1, 3, 3, nT]
  % Both RigidDistances and CorrectRigid accept [nS, 3, 3, nT] or [nS, 9,
  % nT] shapes.
  
  %   for t = 1:nT % ONLY needed for MaxCoilNorms, now done in function call below.
  %     T_Location(:, :, :, t) = ChangeCoordinates(...
  %       reshape(handles.Rigid, [3, 3])', ... % Points are lines
  %       squeeze(T_Location(:, :, :, t))', ... % Points are lines
  %       [], [], [], true)'; % Inverse is true. % [1, 3, 3, nT]
  %   end
  if strcmpi(get(handles.Menu_UseRigid, 'Checked'), 'on') % handles.UseRigid
    % No need to correct for rigid shape here.
    T_LocationDist = squeeze(RigidDistances(T_Location, Reference)); % squeeze([1, 1, nT])
  else
    % Correct the median coil locations to impose rigid body shape.
    T_Location = CorrectRigid(T_Location, handles.Rigid);
    T_LocationDist = squeeze(MaxCoilNorms( bsxfun(@minus, ...
      T_Location, Reference) )); % squeeze([1, 1, nT])
  end
  
  % Distance and FitError quantile of samples per trials at the desired
  % proportion.
  T_P = max(1, min(handles.nS, round((1 - ...
    str2double(get(handles.edit_Proportion, 'String'))/100 ) * handles.nS)));
  T_P2 = max(1, min(handles.nS, round((1 - ...
    str2double(get(handles.edit_Proportion, 'String'))/100 /2 ) * handles.nS)));
  % Only use data from good trials for overall measures.
  nTG = length(handles.GoodTrials);
  P = max(1, min(handles.nS*nTG, round((1 - ...
    str2double(get(handles.edit_Proportion, 'String'))/100 ) * handles.nS*nTG)));
  
  % min and max distances and desired quantile.
  Data = CorrectRigid(handles.Data, handles.Rigid);
  Data = bsxfun(@minus, Data, Reference);
  if strcmpi(get(handles.Menu_UseRigid, 'Checked'), 'on') % handles.UseRigid
    Distances = RigidDistances(handles.Data, Reference); % [nS, 1, nT]
  else
    % Correct coil locations to impose rigid body shape.
    Distances = MaxCoilNorms( Data ); % [nS, 1, nT]
  end
  
  T_Distance_Q = sort(Distances, 1);
  T_MaxDistance = squeeze(T_Distance_Q(handles.nS, :, :)); % squeeze([1, 1, nT])
  T_MinDistance = squeeze(T_Distance_Q(1, :, :));
  T_Distance_Q = squeeze(T_Distance_Q(T_P, :, :)); % nT
  
  MaxAll = max(Distances(:)); % For plot axis limits only.
  DistancesG = Distances(:, :, handles.GoodTrials);
  Distance_Q = sort(DistancesG(:));
  if nTG > 0 % In case of all trials bad.
    set(handles.edit_MaxDistance, 'String', ...
      num2str(Distance_Q(end), handles.NumFormat));
    %   MinDistance = Distance_Q(1); % Expect something close to zero.
    set(handles.edit_Distance_Q, 'String', ...
      num2str(Distance_Q(P), handles.NumFormat));
  else
    Distance_Q = NaN;
    set(handles.edit_MaxDistance, 'String', 'NaN');
    set(handles.edit_Distance_Q, 'String', 'NaN');
  end
  
  % Measures of dispersion: mean distance, "standard deviation", "median
  % absolute deviation".  Still looking for multidimensional dispersion
  % measures...
  %   MeanDistance = mean(MaxMovement(:));
  %   Std = sqrt( 1/(numel(MaxMovement) - 1) * sum(MaxMovement(:).^2) );
  %   MAD = median(MaxMovement(:));
  %   Dispersion = ;
  
  % Upper bound for range of movement.  Max out of the 3 coils of sum of
  % max range in each direction.
  Data = sort(Data, 1);
  T_Range = squeeze(MaxCoilNorms(Data(end, :, :) - Data(1, :, :))); % squeeze([1, 1, nT])
  T_Range_Q = squeeze(MaxCoilNorms(Data(T_P2, :, :) - Data(end+1-T_P2, :, :))); % squeeze([1, 1, nT])
  % Real range can't exceed twice the max distance so enforce that bound.
  T_Range = min(T_Range, 2 * T_MaxDistance);
  Range = reshape(permute( Data(:, :, handles.GoodTrials), ...
    [1, 3, 2]), [handles.nS * nTG, 9]);
  Range = MaxCoilNorms( max(Range, [], 1) - min(Range, [], 1) );
  Range = min(Range, 2 * Distance_Q(end));
  set(handles.edit_Range, 'String', num2str(Range, handles.NumFormat));

  MaxAll = max( MaxAll, max(T_Range(:)) );
  
  
  % ----------------------------------------------------------------------
  % Average coil position fitting error, convert to percent.  Always use
  % max of any coil direction.  Only changes if changing thresholds.
  FitData = max(100 * handles.FitData, [], 2); % [nS, 1, nT]
  T_Fit = squeeze(mean(FitData, 1)); % squeeze([1, 1, nT])
  Fit = mean(T_Fit(handles.GoodTrials));
  
  T_Fit_Q = sort(FitData, 1);
  T_MaxFit = squeeze(T_Fit_Q(end, :, :)); % nT
  T_MinFit = squeeze(T_Fit_Q(1, :, :)); % nT
  T_Fit_Q = squeeze(T_Fit_Q(T_P, :, :)); % nT
  
  FitData = FitData(:, :, handles.GoodTrials);
  %   Fit_Q = sort(FitData(:));
  set(handles.edit_MaxFit, 'String', ...
    num2str(max(FitData(:)), handles.NumFormat)); %Fit_Q(end)
  %   MaxFit = Fit_Q(end);
  %   Fit_Q = Fit_Q(P);
  
  
  % ----------------------------------------------------------------------
  % Mark trials as bad based on thresholds.
  
  Trials = (1:nT)';
  TrialsLoop = [Trials; Trials(end:-1:1)];
  
  % Trials that were already marked bad and not included.
  BadTrials = setdiff(Trials, handles.GoodTrials);
  nTB = length(BadTrials);
  if nTB ~= nT - nTG
    error('Lost trials.');
  end
  
  % Prepare new trial classification but keep old one in case thresholds
  % are changed before saving.
  handles.NewClass = handles.Class;
  % Distance from reference. 'BAD_HeadMotion_Distance'
  % Trial proportion parameter is already taken into account in quantile.
  RejectedTrials = Trials(T_Distance_Q > ...
    str2double(get(handles.edit_Thresh, 'String')) );
  if ~isempty(RejectedTrials)
    handles.NewClass = AddClass(handles.NewClass, 'BAD_HeadMotion_Distance', RejectedTrials);
  end
  
  % Range of motion within trial. 'BAD_HeadMotion_Range'
  RejectedTrials = Trials(T_Range_Q > ...
    str2double(get(handles.edit_Thresh, 'String')) );
  if ~isempty(RejectedTrials)
    handles.NewClass = AddClass(handles.NewClass, 'BAD_HeadMotion_Range', RejectedTrials);
  end
  
  % Fit error. 'BAD_HeadMotion_Fit'
  % Trial proportion parameter is already taken into account in quantile.
  RejectedTrials = Trials(T_Fit_Q > ...
    str2double(get(handles.edit_FitThresh, 'String')) );
  if ~isempty(RejectedTrials)
    handles.NewClass = AddClass(handles.NewClass, 'BAD_HeadMotion_Fit', RejectedTrials);
  end
  
  % Trials that we newly marked as bad, they were included.  Have to loop
  % through all existing classes in case a threshold was made stricter and
  % trials were added.
  NewBadTrials = [];
  for c = 1:numel(handles.NewClass)
    NewBadTrials = [NewBadTrials; handles.NewClass(c).Trials]; %#ok<AGROW>
  end
  NewBadTrials = setdiff(NewBadTrials, BadTrials);
  
  if isempty(NewBadTrials)
    % Nothing to reject, disable button.
    set(handles.button_Reject, 'Enable', 'off');
  else
    % Enable button.
    set(handles.button_Reject, 'Enable', 'on');
  end
    
  % ----------------------------------------------------------------------
  % Draw graphs or update.  Use graph handle to know if already drawn.

  % Set axes ranges.
  Padding = 1.1;
  YMax = Padding * MaxAll;
  if isnan(YMax) % In case of all trials bad.
    YMax = 1;
  end
  YMaxFit = max( min(16, 8 * Fit), Padding * ...
    str2double(get(handles.edit_MaxFit, 'String')) );
  
  % Set proper X limits based on slider values.
  XLim = get(handles.slider_Trials, 'Value') + [0, ...
    nT - get(handles.slider_Trials, 'Max') + 1];
  set(handles.axes_Trials, 'XLim', XLim, 'YLim', [0, YMax], 'Color', 'none');
  set(handles.axes_TrialsFit, 'XLim', XLim, 'YLim', [0, YMaxFit]); %'XTick', [], 'XTickLabel', {''}, ...

  Dark = 153/255;
  Blue = [0, 0, Dark];
  Red = [Dark, 0, 0];
  Green = [0, Dark, 0];
  Grey = [Dark, Dark, Dark]/2;
  Yellow = [1, 1, 0];
  Transp = 0.5;
  Thin = 0.3;
  
  % Highlighting generates a bunch of separate patches, so need to delete
  % and draw again.  Normally only when updating elements, but needed also
  % when restoring the dataset.
  if isfield(handles, 'BadTrPatches')
    delete(handles.BadTrPatches);
    handles.BadTrPatches = [];
  end
  if isfield(handles, 'NewBadTrPatches')
    delete(handles.NewBadTrPatches);
    handles.NewBadTrPatches = [];
  end
  
  if ~isfield(handles, 'Plot_Distance_Q')
    
    % Draw graphs.
    
    set(gcf, 'CurrentAxes', handles.axes_Trials);
    %     axes(handles.axes_Trials); %#ok<*MAXES>
    AxPos = get(handles.axes_Trials, 'Position');
    handles.XLabel = text(AxPos(3)/2, 10 - handles.TextPos, 'Trials', handles.TextOptions{:}, ...
      'Parent', handles.axes_Trials, 'VerticalAlignment', 'bottom');
    handles.YLabel = text(-handles.TextPos, AxPos(4)/2, 'Distance from reference (mm)', ...
      'Parent', handles.axes_Trials, 'Rotation', 90, ...
      'Color', Blue, handles.TextOptions{:}, 'VerticalAlignment', 'top');
    
    % Draw lines first so they are antialiased.
    handles.Plot_LocationDist = plot(Trials, T_LocationDist, ... % 'XData', 'YData',
      'Marker', '.', 'LineStyle', 'none', 'Color', Blue); % 'LineStyle', '-', 'LineWidth', Thick,
    if nTG > 0
      if nT > 1
        handles.Plot_Distance_Q = plot(Trials, T_Distance_Q, ... % 'XData', 'YData',
          'LineStyle', '-', 'LineWidth', Thin, 'Color', Blue); %, ...
      else
        % Plot real distances per sample as a "preview".
        X = 1 + (0:handles.nS-1)'/handles.nS;
        handles.Plot_Distance_Q = plot(X, Distances, ... % 'XData', 'YData',
          'LineStyle', '-', 'LineWidth', Thin, 'Color', Blue); %, ...
      end
    end
    handles.Plot_Thresh = plot([1, nT+1], ... % 'XData'
      [1, 1] * str2double(get(handles.edit_Thresh, 'String')), ... %'YData'
      'LineStyle', ':', 'LineWidth', Thin, 'Color', Blue); %, ...
    handles.Plot_Range_Q = plot(Trials, T_Range_Q, ... % 'XData', 'YData',
      'LineStyle', '-', 'LineWidth', Thin, 'Color', Green); %, ...
    % Draw patches under lines for when transparency is disabled.
    handles.Plot_Range = patch(TrialsLoop, ...
      [T_Range; zeros(nT, 1)], -1 * ones(size(TrialsLoop)), NaN, ... % 'XData', 'YData', 'ZData' (optional), 'CData' (required)
      'LineStyle', '-', 'LineWidth', Thin, ...
      'EdgeColor', Washed(Green), 'FaceColor', Washed(Green), 'FaceAlpha', Transp);
    handles.Plot_DistanceRange = patch(TrialsLoop, ...
      [T_MinDistance; T_MaxDistance(end:-1:1)], -0.5 * ones(size(TrialsLoop)), NaN, ... % 'XData', 'YData', 'ZData' (optional), 'CData' (required)
      'LineStyle', '-', 'LineWidth', Thin, ...
      'EdgeColor', Washed(Blue), 'FaceColor', Washed(Blue), 'FaceAlpha', Transp);
    
    % ___________________________________
    % Second set of axes for fit error.
    set(gcf, 'CurrentAxes', handles.axes_TrialsFit);
    %     axes(handles.axes_TrialsFit);
    handles.YLabelFit = text(AxPos(3) + handles.TextPos, AxPos(4)/2, 'Fit error (%)', ...
      'Parent', handles.axes_TrialsFit, 'Rotation', -90, ...
      'Color', Red, handles.TextOptions{:}, 'VerticalAlignment', 'top');
    
    % Draw lines first so they are antialiased.
    %     handles.Plot_Fit = plot(Trials, T_Fit, ... % 'XData', 'YData',
    %       'XDataSource', 'Trials', 'YDataSource', 'T_Fit', ...
    %       'LineStyle', '-', 'LineWidth', Thick, 'Color', Red); %, ...
    if nTG > 0
      if nT > 1
        handles.Plot_Fit_Q = plot(Trials, T_Fit_Q, ... % 'XData', 'YData',
          'LineStyle', '-', 'LineWidth', Thin, 'Color', Red); %, ...
      else
        % Plot real fit error per sample as a "preview".
        handles.Plot_Fit_Q = plot(X, FitData, ... % 'XData', 'YData',
          'LineStyle', '-', 'LineWidth', Thin, 'Color', Red); %, ...
      end
    end
    handles.Plot_FitThresh = plot([1, nT+1], ... % 'XData'
      [1, 1] * str2double(get(handles.edit_FitThresh, 'String')), ... %'YData'
      'LineStyle', ':', 'LineWidth', Thin, 'Color', Red); %, ...
    % Draw patches under lines for when transparency is disabled.
    handles.Plot_MaxFit = patch(TrialsLoop, ...
      [T_MaxFit; T_MinFit(end:-1:1)], -1 * ones(size(TrialsLoop)), NaN, ... % 'XData', 'YData', 'ZData' (optional), 'CData' (required)
      'LineStyle', '-', 'LineWidth', Thin, ...
      'EdgeColor', Washed(Red), 'FaceColor', Washed(Red), 'FaceAlpha', Transp);
    
  else
    % Only need to update existing plot elements.
    set(handles.Plot_FitThresh, 'XData', [1, nT+1], ...
      'YData', [1, 1] * str2double(get(handles.edit_FitThresh, 'String')) );
    set(handles.Plot_MaxFit, 'XData', TrialsLoop, ...
      'YData', [T_MaxFit; T_MinFit(end:-1:1)], ...
      'ZData', -1 * ones(size(TrialsLoop)) );
    
    set(handles.Plot_LocationDist, 'XData', Trials, 'YData', T_LocationDist );
    if nTG > 0
      if nT > 1
        set(handles.Plot_Distance_Q, 'XData', Trials, 'YData', T_Distance_Q );
        set(handles.Plot_Fit_Q, 'XData', Trials, 'YData', T_Fit_Q);
      else
        % Plot real distances per sample as a "preview".
        X = 1 + (0:handles.nS-1)'/handles.nS;
        set(handles.Plot_Distance_Q, 'XData', X, 'YData', Distances );
        set(handles.Plot_Fit_Q, 'XData', X, 'YData', FitData);
      end
    else
      % Zero existing lines.
      set(handles.Plot_Distance_Q, 'YData', zeros(size(get(handles.Plot_Distance_Q, 'XData'))));
      set(handles.Plot_Fit_Q, 'YData', zeros(size(get(handles.Plot_Fit_Q, 'XData'))));
    end
    set(handles.Plot_Thresh, 'XData', [1, nT+1], ... % 'XData'
      'YData', [1, 1] * str2double(get(handles.edit_Thresh, 'String')) );
    set(handles.Plot_Range_Q, 'XData', Trials, 'YData', T_Range_Q );
    set(handles.Plot_Range, 'XData', TrialsLoop, ...
      'YData', [T_Range; zeros(nT, 1)], ...
      'ZData', -1 * ones(size(TrialsLoop)) );
    set(handles.Plot_DistanceRange, 'XData', TrialsLoop, ...
      'YData', [T_MinDistance; T_MaxDistance(end:-1:1)], ...
      'ZData', -0.5 * ones(size(TrialsLoop)) );
    
  end % draw or update
  
  
  % On second axes so highlighting is under everything.
  set(gcf, 'CurrentAxes', handles.axes_TrialsFit);
  %   axes(handles.axes_TrialsFit);
  
  % Trials that were already marked bad and not included, shade in grey.
  if ~isempty(BadTrials)
    [X, Y] = MakeRectangles(BadTrials, YMaxFit);
    handles.BadTrPatches = patch(X, Y, -2 * ones(size(X)), NaN, ... % 'XData', 'YData', 'ZData' (optional), 'CData' (required)
      'LineStyle', '-', 'LineWidth', Thin, ...
      'EdgeColor', Washed(Grey), 'FaceColor', Washed(Grey), 'FaceAlpha', Transp);
  end
  if ~isempty(NewBadTrials)
    % Highlight new bad trials in yellow.
    [X, Y] = MakeRectangles(NewBadTrials, YMaxFit);
    handles.NewBadTrPatches = patch(X, Y, -2 * ones(size(X)), NaN, ... % 'XData', 'YData', 'ZData' (optional), 'CData' (required)
      'LineStyle', '-', 'LineWidth', Thin, ...
      'EdgeColor', Washed(Yellow), 'FaceColor', Washed(Yellow), 'FaceAlpha', Transp);
  end
  
  % Back to first axes for proper stacking.
  %   set(gcf, 'CurrentAxes', handles.axes_Trials);
  axes(handles.axes_Trials); %#ok<*MAXES>
  
  % ___________________________________
  % Update handles structure to save handles of graph elements and class
  % structure.
  guidata(handles.figure_HeadMotion, handles);
  
  % Reset mouse pointer.
  pause(0.05); % required so that 'BusyAction' 'cancel' works!
  set(handles.figure_HeadMotion, 'Pointer', 'arrow');
end



% ----------------------------------------------------------------------
% Lighten colors.
function Color = Washed(Color)
  Wash = 0.3;
  Color = 1 - (Wash * (1 - Color));
end




% ----------------------------------------------------------------------
% Prepare rectangular regions for plotting shading/highlighting of trials.
function [X, Y] = MakeRectangles(BadTrials, YMax)
  if ~isempty(BadTrials)
    BadTrRangesIx = find(diff(BadTrials) > 1)';
    X = [ BadTrials([1, BadTrRangesIx + 1])' ; ... %- 0.5
      BadTrials([BadTrRangesIx, end])' + 1 ]; %+ 0.5 
    clear BadTrRangesIx
    X = X([1, 1, 2, 2], :);
    Y = YMax * [0; 1; 1; 0] * ones([1, size(X, 2)]);
  else
    X = [];
    Y = [];
  end
end





