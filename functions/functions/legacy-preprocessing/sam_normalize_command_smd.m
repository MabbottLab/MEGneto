function sam_normalize_command_smd(mrilist, samlist, num, res)

% sam_normalize_command_smd(mrilist, samlist, num, res)
% mrilist is path to single *.MRI file
% samlist is path to single *.IMG file converted from *.SVL derived from
%         headmodel
% num is number of subjects
% res is resolution in mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAM_normalize 
% gui to select multiple analyze format SAM images to normalize them using
% the T1 template with brain masking in SPM2
%
% -- assumes that co-registered .mri and .svl files have been converted
% to analyze format using the svl2spm program
%
% D. Cheyne
% June, 2005
%
%	revisions:
%	July 5, 2005	-- fix to force defaults.analyze.flip = 0 so that
%			   RAS format analyze files (from svl2spm) are not
%			   flipped during normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

versionNo = '1.1';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% configure parameters for normalization

spm_defaults;

defaults.analyze.flip = 0;	% input files from svl2spm are RAS -- SPM2 

% Spatial Normalisation defaults
%=======================================================================
% defaults.normalise.estimate.smosrc  = 8;
% defaults.normalise.estimate.smoref  = 0;
% defaults.normalise.estimate.regtype = 'mni';
% defaults.normalise.estimate.weight  = '';
% defaults.normalise.estimate.cutoff  = 25;
% defaults.normalise.estimate.nits    = 16;
% defaults.normalise.estimate.reg     = 1;
% defaults.normalise.estimate.wtsrc   = 0;
% defaults.normalise.write.preserve   = 0;
% defaults.normalise.write.bb         = [[-78 -112 -50];[78 76 85]];
% defaults.normalise.write.vox        = [2 2 2];
% defaults.normalise.write.interp     = 1;
% defaults.normalise.write.wrap       = [0 0 0];

%% get current default values for normalization and write to output
defs = defaults.normalise;

%% override defaults so we know explicitly what is being used
flags = defs.estimate;

flags.nits = 12;         % 12 is enough ?
%flags.reg = 0.1;       

estmation_params = flags;

write_flags = defs.write;

% write_flags.bb = [[-78 -112 -50];[78 76 85]];
% write_flags.wrap = [0 0 0];
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



spm('FnBanner',mfilename, versionNo);

%spm_figure('Clear','Interactive');
[Finter,~,CmdLine] = spm('FnUIsetup','SAM Image normalization');

%%
%% get common parameters
%%
numSubjects = num; %spm_input('# Subjects to normalize'); 
if ( numSubjects < 1 )
    sprintf('exiting...\n')
   return;
end

%res = spm_input('Resolution (in mm) of normalized images'); 
imageResolution = [res res res];

write_flags.vox = imageResolution;

write_params = write_flags;

%%
%% get list of names
%%
% for i=1:numSubjects
% 
%     % select structural
%     label = sprintf('Subject # %d: Select MRI file (.img) to use for normalization', i);
%     mriList(i).name = spm_get(Inf,'.img',label);  
%     
% 	if isempty(mriList(i).name)
%         return; 
%     end
% 
%     % select all SAM files at once...
%    label = sprintf('Subject # %d: Select SAM(.img) file(s) to  normalize', i);
%    samList(i).names = spm_get(Inf,'.img',label);  
% 
% end;
 
%% 
%% loop through all subjects and images
%%
mriList = mrilist;
samList = samlist;

for i=1:numSubjects
    
    disp(['Processing subject #' num2str(i)]);

  
    imageToGetNormalizationFrom = mriList;
    
    samFiles = samList;
    
    numImages=size(samFiles,1);
   
    for j=1:numImages
        
        
        imageToNormalize = deblank(samFiles(j,:));
 
    
        %%
        %% need to compute normalization parameters for this MRI once only 
        %%
    
        if ( j == 1 )
        
            templateFile = fullfile(spm('Dir'),'templates/T1.mnc');
            maskFile = fullfile(spm('Dir'),'apriori/brainmask.mnc');
            objMaskFile = '';
            
            [mriPath,mriRoot,mriExt]=fileparts(imageToGetNormalizationFrom);        
            
            snMatFile= sprintf('%s/%s_%gmm_sn3d.mat', mriPath,mriRoot,res);   % mat file to store params
      
            disp(['']);
            spm('FigName',['Finding Norm params for subject ' num2str(i)],Finter,CmdLine);
            disp('------------------------------------');
            disp(['Calculating normalisation parameters from ' imageToGetNormalizationFrom ]);

            flags.graphics = 0;     % spm_normalise crashes on display code, so turn display graphics off
            spm_normalise(templateFile, imageToGetNormalizationFrom, snMatFile, maskFile, objMaskFile, flags);
            disp(['Saving normalisation parameters in ' snMatFile ]);       
        end
     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % apply normalization to functional images
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        disp(['Applying normalisation to ' imageToNormalize ]);
        spm_write_sn( imageToNormalize, snMatFile, write_flags);
    
        disp('------------------------------------');
        
    end  % next subject

end;

disp('all finished!!!');
spm('FigName','normalization all finished!',Finter,CmdLine);

