function [hF,hB,fig] = imoverlay(B,F,climF,climB,cmap,alpha,haxes,nslices,range,title)
% IMOVERLAY(B,F) displays the image F transparently over the image B.
%    If the image sizes are unequal, image F will be scaled to the aspect
%    ratio of B.
% 
%    [hF,hB] = imoverlay(B,F,[low,high]) limits the displayed range of data
%    values in F. These values map to the full range of values in the
%    current colormap.
% 
%    [hF,hB] = imoverlay(B,F,[],[low,high]) limits the displayed range of
%    data values in B.
% 
%    [hF,hB] = imoverlay(B,F,[],[],map) applies the colormap to the figure.
%    This can be an array of color values or a preset MATLAB colormaps
%    (e.g. 'jet' or 'hot').
% 
%    [hF,hB] = imoverlay(B,F,[],[],[],alpha) sets the transparency level to
%    alpha with the range 0.0 <= alpha <= 1.0, where 0.0 is fully
%    transparent and 1.0 is fully opaque.
% 
%    [hF,hB] = imoverlay(B,F,[],[],[],[],ha) displays the overlay in the
%    axes with handle ha.
% 
%    [hF,hB] = imoverlay(B,F,[],[],[],[],[],slice_range) displays slices
%    specified by slice_range. Default is 1:80.
%
%    [hF,hB] = imoverlay(B,F,[],[],[],[],[],[],nslice) determines
%    number of slice to display. Default every slice in slice_range.
%
%    [hF,hB] = imoverlay(...) returns the handles to the front and back
%    images.
%
%
% Author: Matthew Smith / University of Wisconsin / Department of Radiology
% Date created:  February 6, 2013
% Last modified: Jan 2, 2015
% EDITED BY DIANA.
%
%  
%  Examples:
%     
%     % Overlay one image transparently onto another
%     imB = phantom(256);                       % Background image
%     imF = rgb2gray(imread('ngc6543a.jpg'));   % Foreground image
%     [hf,hb] = imoverlay(imB,imF,[40,180],[0,0.6],'jet',0.6);
%     colormap('parula'); % figure colormap still applies
%
%
%     % Use the interface for flexibility
%     imoverlay_tool;
%
% 
% See also IMOVERLAY_TOOL, IMAGESC, HOLD, CAXIS.



ALPHADEFAULT = 0.4; % Default transparency value
CMAPDEFAULT = 'parula';

if nargin == 0,
    try
        imoverlay_tool;
        return;
    catch
        errordlg('Cannot find imoverlay_tool.', 'Error');
    end
end


% Check image sizes
if ~isequal(size(B),size(F))
    fprintf('Warning! Image sizes unequal. Undesired scaling may occur.\n');
end

% Check arguments
if nargin < 10 || isempty(title)
    title = 'BSR Thresholded Salience Strength';
end

if nargin < 9 || isempty(range)
    range = 1:80;
end

if nargin < 8 || isempty(nslices)
    nslices = length (range);
end

if nargin < 7 || isempty(haxes)
    haxes = [];          
end

if nargin < 6 || isempty(alpha)
    alpha = ALPHADEFAULT;
end

if nargin < 5 || isempty(cmap)
    cmap = CMAPDEFAULT;    
end

if nargin < 4 || isempty(climB)
    climB = [min(B(:)), max(B(:))];
end

if nargin < 3 || isempty(climF)
    climF = [min(F(:)), max(F(:))]; 
end

if abs(alpha) > 1
    error('Alpha must be between 0.0 and 1.0!');
end

% determine which slices to use
if nslices == length(range)
    index = range;
else
    index = round(linspace(range(1),range(end),nslices));
end

% Create a figure unless axes is provided
if isempty(haxes) || ~ishandle(haxes)
%     f=figure('Visible','off',...
%         'Units','pixels','Renderer','opengl');
%     pos = get(f,'Position');
%     set(f,'Position',[pos(1),pos(2),size(B,2),size(B,1)]);
%     haxes = axes;
%     set(haxes,'Position',[0,0,1,1]);
%     movegui(f,'center');
    if range(1) == 1 && range(end) == 80 && nslices == 80 ;
        fig = figure('units','pixels','Position',[200 100 1090 728],'Name',title); 
        [ha,pos] = tight_subplot(8,10,0,0,0);
    else % determine dimension of subplots needed (only works for ranges divisible by integers)
        mod_array=zeros(1,12);
        for i=1:12
            mod_array(i)=mod(nslices,i);
        end
        if nnz(mod_array==0) > 1 
            [~,j1] = find(mod_array==0);
            j2=nslices./j1;
            [~,centre_value] = min(abs(j2-j1));
        else
            error('Range of slices must be divisible by whole numbers.');
        end
        fig = figure('units','pixels','Position',[200 100 109*j2(centre_value) 91*j1(centre_value)],...
            'Name',title);
        [ha,pos] = tight_subplot(j1(centre_value),j2(centre_value),0,0,0);
    end

    % Create colormap
    cmapSize = 100; % default size of 60 shows visible discretization
    if ischar(cmap)
        
        try
            cmap = eval([cmap '(' num2str(cmapSize) ');']);
        catch
            fprintf('Colormap ''%s'' is not supported. Using ''jet''.\n',cmap);
            cmap = jet(cmapSize);
        end
    end
    colormap(cmap);
end



for ss = 1:nslices;
    
    new_B = B(:,:,index(ss));    
    new_F = F(:,:,index(ss));

    % To have a grayscale background, replicate image to 3-channels
    new_B = repmat(mat2gray(double(new_B),double(climB)),[1,1,3]);

    % Display the back image
    axes(ha(ss));
    hB = imagesc(new_B,climB);  axis image off;
    set(gca,'Position',pos{ss});

    % Add the front image on top of the back image
    hold on;
    axes(ha(ss));
    hF = imagesc(new_F,climF);
    set(gca,'Position',pos{ss});


    % Make the foreground image transparent
    alphadata = alpha.*(new_F >= climF(1));
    set(hF,'AlphaData',alphadata);    


end

end