function [] = sam_svl2spm_command_smd(mrilist, svllist, num)

% sam_svl2spm_command(mrilist, svllist, num)
% mrilist(ii).names contains full path to mri
% svllist(ii).names contains full path of each svl file to be normalized
% If there is more than one svl file per mri they should be organized into
% an mXn matrix where m = num of files/mri and n is length of longest
% filename
% num is now many subjects with mri

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SAM_svl2spm.m
%   
%  gui for conversion of CTF SAM (.svl) and corresponding MRI (.mri) files
%  to Analyze format for SPM. 
%  the SAM image is saved as axial (neurological) Analyze format
%  (.img and .hdr).
%  the MRI image is resliced to fit the bounding box of the SAM image and
%  saved in the same format (mriFilename_resl.img and
%  mriFilename_resl.hdr).
%
%  Dependencies:
%                  - SPM2 toolbox 
%                  (http://www.fil.ion.ucl.ac.uk/spm/software/spm2/)
%                  - MRI/EEG - toolbox (http://eeg.sourceforge.net/)
%                      - ctf_read_mri
%                      - avw_hdr_make
%                      - avw_write (slightly modified not to define SPM
%                      origin!!!)
%                  - SPM2 SAM toolbox
%                      - ctf_svl2spm
%                      - ctf_svl_mriReslice
%
%  Versions:
%                  - Version 1.0: used svl2spm Version 1.3 written by 
%                  D. Cheyne, May 2005 (C++ code); written by D. Cheyne.
%                  - Version 1.1: modified to use Matlab functions
%                  (ctf_svl2spm and ctf_svl_mriReslice) adaptations of
%                  svl2spm; written February 2006
%                  (Andreea Bostan, Hospital for Sick Children) 
%
%  
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

versionNo = '1.1';
spm('FnBanner',mfilename, versionNo);

%%%
% SPM interface
%%%
[Finter,~,CmdLine] = spm('FnUIsetup','SAM to SPM conversion');

%%%
% Input number of files to convert 
%%%
numSubjects = num;%spm_input('# Subjects with files to convert'); 
if ( numSubjects < 1 )
    sprintf('exiting...\n')
   return;
end

%%%
% Type of interpolation
%%%
interp_method = 2; %spm_input('Interpolation:',2,'b','none|nearest|linear|cubic', [0,1,2,3],3)

%%%
% Get list of files
%
% for i=1:numSubjects
% 
%     %Get structural CTF MRI (.mri) file
%     label = sprintf('Subject# %d: Select MRI (.mri)', i);
%     mriList(i).name = spm_get(Inf,'.mri',label);  
%     
% 	if isempty(mriList(i).name)
%         return; 
%     end % if
% 
%     %Select SAM (.svl) files all at once
%    label = sprintf('Subject# %d: Select SAM (.svl) file(s)', i);
%    svlList(i).names = spm_get(Inf,'.svl',label);  
% 
% end; % numSubjects

mriList = mrilist;
svlList = svllist;

%%%
%Loop through all subjects and convert images
%%%
for i=1:numSubjects
    
    disp(['Processing subject #' num2str(i)]);

    mriFile = mriList;
    
    svlFiles = svlList;
    
    numImages = size(svlFiles,1);
   
    for j=1:numImages
        
        
        svlFile = deblank(svlFiles(j,:));
        disp(['Converting ' svlFile ]);
        
        [~, xVoxels, yVoxels, zVoxels, svlResolution, xmin, xmax, ymin, ymax, zmin, zmax] = ctf_svl2spm(svlFile);

        %unix(['svl2spm' ' -m ' mriImage ' -s ' imageToConvert]); % svl2spm
   
        disp('------------------------------------');
        
    end % numImages
    
    disp(['Converting ' mriFile ]);
    ctf_svl_mriReslice(mriFile, xVoxels, yVoxels, zVoxels, svlResolution, xmin, xmax, ymin, ymax, zmin, zmax, interp_method);

end % numSubjects

disp('All Conversions are Complete !!!');
spm('FigName','DONE',Finter,CmdLine);

return
