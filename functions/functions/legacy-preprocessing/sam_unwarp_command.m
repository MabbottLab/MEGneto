function coords = sam_unwarp_command(sn3d, img, coordsFile, talUnits)
%
% coords = sam_unwarp_command(sn3d, img, coordsFile, talUnits)
%
% [sn3d] is the sn3d.mat warping file
% [img] is an unwarped ANALYZE (*.img) file 
% [coordsFile] is mx3 matrix of mni or tlc coords. m=how many points.
% [talUnits]=1 for tlrc, 0 for mni
% [coords]=tlrc or mni output coords, structure array with field "subject"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Inputs:  - sn3d.mat - mat file created during spatial normalization
%            - .img     - original (unwarped) SAM image in Analyze format
%            - coordinate(s) in Talairach or MNI space
%             
%   
%   Notes: 
%   - Uses SPM2 plus get_orig_coord2() from John Ashburner and mni2tal() 
%   - Images must have been normalized with version 1.2 of svl2spm()
%
%   written by D. Cheyne, July, 2005
%
%   revisions
%   
%      1.1  - removed text output that assumed the large SPM window was open
%
%      1.2  - put in fix for tal2mni bug for 3 x 3 input matrices - needs
%      testing!!!!




spm_defaults;			% load SPM defaults
defaults.analyze.flip = 0;	% force default to don't left-right flip images
versionNo = 1.2;
versionStr = sprintf('SAM_unwarp(version %.1f)', versionNo);

[m n] = size(img);

for zz = 1:m

    imgfile = img(zz,:);
    matfile = sn3d(zz,:);
    coord = coordsFile;

    V = spm_vol(imgfile);

    s = V.private.hdr.hist.aux_file;
    origin = sscanf(s, '%*s %f %f %f');
    res = V.mat(1,1);


    A=size(coordsFile);
    numPts = A(1);

    if (talUnits == 1)

        % fix for 3 x 3 bug in tal2mni -- this code _should_ work but didn't test
        % it yet
        for i=1:numPts
            tempcoord = coord(i,:);
            coord(i,:) = tal2mni(tempcoord);
        end

    end


    % get voxels

    orig_coord = get_orig_coord2(coord, matfile, imgfile);

    % subtract one to get to svl indices

    orig_coord = round((orig_coord-1)*100000)/100000;

    % convert to CTF coords

    x = ( orig_coord(:,2) - origin(2) ) .* res;  % swap x and y
    y = ( origin(1) - orig_coord(:,1) ) .* res;
    z = ( orig_coord(:,3) - origin(3) ) .* res;
    ctf_coord_in_cm = [x y z] ./ 10.0;


    disp(' ');
    t1 = sprintf('Unwarped points for %s:', imgfile);
    disp(t1);
    t1 = sprintf('MNI coordinates(mm)\t  Original CTF coordinates (cm)');
    disp(t1);
    disp('-------------------------------------------------------');

    for i=1:numPts
        t1 = sprintf('%5.1f %5.1f %5.1f ==>\t%5.1f %5.1f %5.1f', ...
            coord(i,1), coord(i,2), coord(i,3), ...
            ctf_coord_in_cm(i,1), ctf_coord_in_cm(i,2), ctf_coord_in_cm(i,3) );
        coords(zz).subject(i, :) =[ctf_coord_in_cm(i,1) ctf_coord_in_cm(i,2) ctf_coord_in_cm(i,3)];
        disp(t1);
    end
end