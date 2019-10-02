%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  [mriImg] = ctf_svl_mriReslice(mriFile, xVoxels, yVoxels, zVoxels, svlResolution, xmin, xmax, ymin, ymax, zmin, zmax, interp_method)
%%
%%  reads and reslices CTF MRI (.mri) file to match a corresponding .svl image
%%  and saves the MRI in Analyze format. 
%%  the MRI image is resliced to fit the bounding box of the SAM (.svl)
%%  image. 
%%  
%%
%%  Inputs:
%%                  - mriFile: .mri filename (string): 'filename.mri'.
%%                  - xVoxels, yVoxels, zVoxels: .svl file dimensions 
%%                  (e.g., from ctf_svl2spm)
%%                      - xVoxels: [coronal]: posterior - anterior.
%%                      - yVoxels: [sagittal]:  left - right.
%%                      - zVoxels: [axial]: inferior - superior.
%%                  - svlResolution: .svl file resolution 
%%                  (e.g. from ctf_svl2spm)
%%                  - xmin, xmax, ymin, ymax, zmin, zmax: X, Y, Z start &
%%                  end (mm) (e.g., from ctf_svl2spm).
%%                  - interp_method: interpolation method to be used for
%%                  resamping of MRI data from original .mri file
%%                      - 'none': no interpolation is used, coefficients
%%                      are rounded.
%%                      - 'nearest': nearest neighbour interpolation
%%                      - 'linear': linear interpolation - recommended
%%                      - 'cubic': cubic interpolation
%%
%%  Output:
%%                  - mriImg: matlab struct with .mri file info and data
%%
%%  Contributions:
%%                  - reading of the .mri file: ctf_read_mri, version 1.5, 
%%                  ritten by D. Weber.
%%                  (http://eeg.sourceforge.net/)
%%                  - reslicing of .mri file: svl2spm, version 1.3, written
%%                  by D. Cheyne (C++ code).
%%
%%  Dependencies:
%%                  - SPM2 toolbox 
%%                  (http://www.fil.ion.ucl.ac.uk/spm/software/spm2/)
%%                  - MRI/EEG - toolbox (http://eeg.sourceforge.net/)
%%                      - ctf_read_mri
%%                      - avw_hdr_make
%%                      - avw_write (slightly modified not to define SPM
%%                      origin!!!)
%%
%%  Versions:
%%                  - Version 1.0: written February 2006 
%%                  (Andreea Bostan, Hospital for Sick Children)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mriImg] = ctf_svl_mriReslice(mriFile, xVoxels, yVoxels, zVoxels, svlResolution, xmin, xmax, ymin, ymax, zmin, zmax, interp_method)

version = 1.0;

fprintf('ctf_svl_mriReslice version %.1f\n\n', version);

mri = ctf_read_mri(mriFile); % read .mri file (=> matlab struct)

%% rotation matrix
fprintf('MRI to head coordinates transformation matrix:\n');
RMat = mri.hdr.transformMatrixHead2MRI(1:3, 1:3);
RMat'

%% origin
fprintf('Origin:\n');
headOrigin = [ mri.hdr.headOrigin_sagittal mri.hdr.headOrigin_coronal mri.hdr.headOrigin_axial ]

%% Voxel Size
if ((mri.hdr.mmPerPixel_sagittal == mri.hdr.mmPerPixel_coronal) & (mri.hdr.mmPerPixel_sagittal == mri.hdr.mmPerPixel_axial))
    mmPerVox = mri.hdr.mmPerPixel_sagittal;
else
    fprintf('ERROR: Voxels dimensions in MRI file do not match!'); return
end


%% all values in mm
%% we need to reslice mri.img in analyze orientation
%% make 256 x 256 axial slices from inferior to superior at 1mm resolution

sagRes = (yVoxels * svlResolution) / 256.0;
corRes = (xVoxels * svlResolution) / 256.0;

numSlices = round(zVoxels * svlResolution)

%% build head location matrix (3x(256*256*numSlices) array)
fprintf('Calculating head location matrix ...\n');
tic
y = ymax - (0:255) .* sagRes;
x = xmin + (0:255) .* corRes; % X values
xy = [reshape(repmat(x, 256, 1), 1, 256*256); repmat(y, 1, 256)]; % tile and reshape XY values
clear x y;
z = zmin + (0:(numSlices-1)); % Z values

HeadLoc = [repmat(xy, 1, numSlices); reshape(repmat(z, 256*256, 1), 1, 256*256*numSlices)]; 
clear z;
clear xy;
toc

switch interp_method
    case 0
        interp_method = 'none';
    case 1
        interp_method = 'nearest';
    case 2
        interp_method = 'linear';
    case 3
        interp_method = 'cubic';
end

%% rotate and scale
fprintf('Converting into MR space ...\n');
tic
MriLoc = (1/mmPerVox) * (RMat * HeadLoc); % same size as HeadLoc
if (interp_method == 0)
    MriVox = int32(MriLoc + repmat(headOrigin', 1, 256*256*numSlices)); % round
else
    MriVox = MriLoc + repmat(headOrigin', 1, 256*256*numSlices); % no rounding
end

clear MriLoc;
toc

%% create final image
fprintf('Creating image ...\n');
tic

if (interp_method == 0)
    
    %% create a not-out-of-bonds mask
    noob = (MriVox > 0 & MriVox <= 256);
    noob = (noob(1,:) & noob(2,:) & noob(3,:));

    Img = zeros(256, 256, numSlices);
    data_ind = sub2ind(size(mri.img), MriVox(1,noob), MriVox(2,noob), MriVox(3,noob));
    Img(noob) = mri.img(data_ind);
    clear data;
    clear data_ind;
    
else
           
    switch interp_method
        case 1
            interp_method = 'nearest';
        case 2
            interp_method = 'linear';
        case 3
            interp_method = 'cubic';
    end
    
    Img = interp3(mri.img, reshape(MriVox(2, :), 256, 256, numSlices), reshape(MriVox(1, :), 256, 256, numSlices), reshape(MriVox(3, :), 256, 256, numSlices), interp_method,0);
    
end

if (mri.hdr.dataSize == 2)
    Img = (1/256) * Img;
end; % if original .mri was 16 bit, scale to 8 bit
toc

%% create default Analyze header and update relevant fields
mriImg = avw_hdr_make;

mriImg.fileprefix = strrep(mriFile, '.mri', '_resl');

%% dimensions and resolution
mriImg.hdr.dime.dim(2) = 256;
mriImg.hdr.dime.dim(3) = 256;
mriImg.hdr.dime.dim(4) = numSlices;
mriImg.hdr.dime.pixdim(2) = sagRes;
mriImg.hdr.dime.pixdim(3) = corRes;
mriImg.hdr.dime.pixdim(4) = 1;
mriImg.hdr.dime.pixdim(5:8) = [0 0 0 0];

mriImg.img = Img;


%% datatype
mriImg.hdr.dime.bitpix = 8;
mriImg.hdr.dime.datatype = 2;


%% write file
avw_write(mriImg, mriImg.fileprefix);

return

