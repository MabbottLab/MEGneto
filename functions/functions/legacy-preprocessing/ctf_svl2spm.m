%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  [svlImg, xVoxels, yVoxels, zVoxels, svlResolution, xmin, xmax, ymin, ymax, zmin, zmax] = ctf_svl2spm(svlFile)
%%
%%  reads and converts CTF SAM (.svl) files to Analyze (.img) format for SPM
%%  and mri3dX. 
%%  The SAM image is saved as axial (neurological) format Analyze files
%%  (.img and .hdr).
%%
%%  Input:
%%                  - svlFile: .svl filename (string): 'filename.svl'.
%%
%%  Outputs:
%%                  - svlImg: matlab struct with .svl file info and data
%%                  - xVoxels, yVoxels, zVoxels: .svl file dimesions
%%                     - xVoxels: [coronal]: posterior - anterior.
%%                     - yVoxels: [saggital]:  left - right.
%%                     - zVoxels: [axial]: inferior - superior. 
%%                  - xmin, xmax, ymin, ymax, zmin, zmax: X, Y, Z start &
%%                  end (mm)
%%
%%  
%%  Contributions:
%%                  - reading of the .svl file: ctf_read_svl, version 1.2,
%%                  written by S. Dalal and D. Weber. 
%%                  (http://eeg.sourceforge.net/)
%%                  - conversion to Analyze format: svl2spm, version 1.3,
%%                  written by D. Cheyne (C++ code).
%%  
%%  Dependencies:
%%                  - SPM2 toolbox 
%%                  (http://www.fil.ion.ucl.ac.uk/spm/software/spm2/)
%%                  - MRI/EEG - toolbox (http://eeg.sourceforge.net/)
%%                      - avw_hdr_make
%%                      - avw_write (slightly modified not to define SPM
%%                      origin!!!)
%%
%%  Versions:
%%                  - Version 1.0: written February 2006                
%%                  (Andreea Bostan, Hospital for Sick Children)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [svlImg, xVoxels, yVoxels, zVoxels, svlResolution, xmin, xmax, ymin, ymax, zmin, zmax] = ctf_svl2spm(svlFile)

version = 1.0;

fprintf('ctf_svl2spm version %.1f\n\n', version);

fid = fopen(svlFile, 'r', 'b', 'latin1');%latin1 was added as a matlab upgrade (to 7.4) required it

%%  read .svl header
identity = transpose(fread(fid,8,'*char'));
if(~strcmp(identity,'SAMIMAGE'))
    error('This doesn''t look like a SAM IMAGE file.');
end % if SAM image
vers = fread(fid,1,'int32'); % SAM file version
setname = fread(fid,256,'*char');
numchans = fread(fid,1,'int32');
numweights = fread(fid,1,'int32');
if(numweights ~= 0)
    warning('... numweights ~= 0');
end

padbytes1 = fread(fid,1,'int32');

XStart = fread(fid,1,'double');
XEnd = fread(fid,1,'double');
YStart = fread(fid,1,'double');
YEnd = fread(fid,1,'double');
ZStart = fread(fid,1,'double');
ZEnd = fread(fid,1,'double');
StepSize = fread(fid,1,'double');

hpFreq = fread(fid,1,'double');
lpFreq = fread(fid,1,'double');
bwFreq = fread(fid,1,'double');
meanNoise = fread(fid,1,'double');

MRIname = transpose(fread(fid,256,'*char'));
nasion = fread(fid,3,'int32');
rightPA = fread(fid,3,'int32');
leftPA = fread(fid,3,'int32');

SAMtype = fread(fid,1,'int32');
SAMunit = fread(fid,1,'int32');
 
padbytes2 = fread(fid,1,'int32');

if ( vers > 1 )
    nasion_meg = fread(fid,3,'double');
    rightPA_meg = fread(fid,3,'double');
    leftPA_meg = fread(fid,3,'double');
    SAMunitname = fread(fid,32,'*char');
end % version 2 has extra fields

SAMimage = fread(fid,inf,'double'); % 1-d array of voxel values
        
fclose(fid);

%%  .svl file is a stack of coronal slices
xVoxels = size(XStart:StepSize:XEnd,2); % posterior -> anterior (coronal) 
yVoxels = size(YStart:StepSize:YEnd,2); % right -> left (saggital)
zVoxels = size(ZStart:StepSize:ZEnd,2); % bottom -> top (axial)

numImageVoxels = xVoxels * yVoxels * zVoxels;
svlResolution = StepSize * 1000.0;
xmin = XStart * 1000.0;
ymin = YStart * 1000.0;
zmin = ZStart * 1000.0;
xmax = XEnd * 1000.0;
ymax = YEnd * 1000.0;
zmax = ZEnd * 1000.0;

fprintf('svl dimesions %d (cor) x %d (sag) x %d (axi) (res = %g mm)\n', xVoxels, yVoxels, zVoxels, svlResolution);

%% create Analyze header and update relevant fields
svlImg = avw_hdr_make; % creates default header structure

svlImg.fileprefix = strrep(svlFile, '.svl', ''); % file prefix

%% dimensions and resolution
svlImg.hdr.dime.dim(2) = yVoxels; 
svlImg.hdr.dime.dim(3) = xVoxels;
svlImg.hdr.dime.dim(4) = zVoxels;
svlImg.hdr.dime.pixdim(2) = svlResolution;
svlImg.hdr.dime.pixdim(3) = svlResolution;
svlImg.hdr.dime.pixdim(4) = svlResolution;
svlImg.hdr.dime.pixdim(5:8) = [0 0 0 0];

%%  format (Analyze) SAM image for writing
%%  .svl array is [posterior -> anterior (coronal) slices] [right -> left ]
%%  [bottom -> top]
%%  Analyze default axial orientation is actually radiological, but we want
%%  neurological orientation to be compatible to how mri3dX reads in
%%  normalized images [left -> right]
%%  we want [axial slices (slowest)] [posterior -> anterior] [left -> right
%%  (fastest)]
Img = reshape(SAMimage, zVoxels, yVoxels, xVoxels); % reshape 1-d array to 3-d
Img = permute(Img, [2 3 1]); % Analyze format
Img = flipdim(Img, 1); % left -> right
svlImg.img = Img;

%%  get min max values of voxel values; these are re-calculated by avw_write
peakNegValue = min(SAMimage);
peakPosValue = max(SAMimage);

svlImg.hdr.dime.glmax = peakPosValue;
svlImg.hdr.dime.glmax = peakNegValue;

%%  datatype ('DOUBLE')
svlImg.hdr.dime.datatype = 64;
svlImg.hdr.dime.bitpix = 64;

%%  get Analyze origin in voxels -- this is relative to RAS origin!
xOriginVox = ymax / svlResolution;
yOriginVox = -xmin / svlResolution;
zOriginVox = -zmin / svlResolution;

fprintf('Analyze origin (in voxels) %.1f %.1f %.1f \n', xOriginVox, yOriginVox, zOriginVox);

origStr = char(zeros(1, 24));
Str = sprintf('ORIGIN %.1f %.1f %.1f', xOriginVox, yOriginVox, zOriginVox);
origStr(1:size(Str, 2)) = Str;

%% store origin in .hdr.hist.aux_file for SAM_unwarp
svlImg.hdr.hist.aux_file = origStr;

%%  write file
avw_write(svlImg, svlImg.fileprefix);

return

