function [hF,hB] = imoverlay_timeseries(B,F,T,cmap,alpha,range)

ALPHADEFAULT = 0.4; % Default transparency value
CMAPDEFAULT = 'parula';


% Check image sizes
if ~isequal(size(B),size(F))
    fprintf('Warning! Image sizes unequal. Undesired scaling may occur.\n');
end

% Check arguments
if nargin < 6 || isempty(range)
    range = 1:80;
end

if nargin < 5 || isempty(alpha)
    alpha = ALPHADEFAULT;
end

if nargin < 4 || isempty(cmap)
    cmap = CMAPDEFAULT;    
end

clim = [abs(min(F(:))),abs(max(F(:)))];
climB = [min(B(:)), max(B(:))];
climF = [-1*(max(clim)) max(clim)];

if abs(alpha) > 1
    error('Alpha must be between 0.0 and 1.0!');
end


% Create a figure unless axes is provided
%     f=figure('Visible','off',...
%         'Units','pixels','Renderer','opengl');
%     pos = get(f,'Position');
%     set(f,'Position',[pos(1),pos(2),size(B,2),size(B,1)]);
%     haxes = axes;
%     set(haxes,'Position',[0,0,1,1]);
%     movegui(f,'center');
if range(1) == 1 && range(end) == 80 ;
    size_range=length(range);
    f = figure('Visible','off','units','pixels','Position',[200 100 1090 728]);
    [ha,pos] = tight_subplot(8,10,[.01 .03],[.2 .01],[.01 .01]);
else % determine dimension of subplots needed (only works for ranges divisible by integers)
    size_range=length(range);
    mod_array=zeros(1,12);
    for i=1:12
        mod_array(i)=mod(size_range,i);
    end
    if nnz(mod_array==0) > 1 
        [~,j1] = find(mod_array==0);
        j2=size_range./j1;
        [~,centre_value] = min(abs(j2-j1));
    else
        error('Range of slices must be divisible by whole numbers.');
    end
    f = figure('Visible','off','units','pixels','Position',[200 100 109*j2(centre_value) 91*j1(centre_value)]);
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


% % initial display
% time = 1;
% atlas_map=zeros(91,109,91);
% for ss=1:size(T,1) 
%     atlas_map(F == ss) = T(ss,time,1);
% end
% 
% for ss = 1:size_range;
%     
%     nwe_B = B(:,:,range(ss));    
%     new_F = atlas_map(:,:,range(ss));
% 
%     % To have a grayscale background, replicate image to 3-channels
%     nwe_B = repmat(mat2gray(double(nwe_B),double(climB)),[1,1,3]);
% 
%     % Display the back image
%     axes(ha(ss));
%     hB = imagesc(nwe_B,climB);  axis image off;
%     set(gca,'Position',pos{ss});
% 
%     % Add the front image on top of the back image
%     hold on;
%     axes(ha(ss));
%     hF = imagesc(new_F,climF);
%     set(gca,'Position',pos{ss});
% 
% 
%     % Make the foreground image transparent
%     alphadata = alpha.*(new_F >= climF(1));
%     set(hF,'AlphaData',alphadata);
% 
% end
f = figure('Visible','off');%,'Units','inches','Position',[2 2 13.5 9]);
% create slider to control time
sld=uicontrol(f,'Style', 'slider',...
    'BackgroundColor', 'c','Min',1,'Max',100,'Value',60) %,...
%     'Units','inches','Position', [3 3 3 0.5],...
%     'Callback', @surfzlim);
% add text to slider
txt= uicontrol(f,'Style','text','String','Time (s)');

%     'Units','inches','Position', [3 3 3 0.5],...
% make figure visible after adding all components
set(f,'Visible','on');

 % Create a figure and axes
    f = figure('Visible','off');
    ax = axes('Units','pixels');
    surf(peaks)
    
    % Create pop-up menu
    popup = uicontrol('Style', 'popup',...
           'String', {'parula','jet','hsv','hot','cool','gray'},...
           'Position', [20 340 100 50],...
           'Callback', @setmap);    
    
   % Create push button
    btn = uicontrol('Style', 'pushbutton', 'String', 'Clear',...
        'Position', [20 20 50 20],...
        'Callback', 'cla');       

   % Create slider
    sld = uicontrol('Style', 'slider',...
        'Min',1,'Max',50,'Value',41,...
        'Position', [400 20 120 20],...
        'Callback', @surfzlim); 
					
    % Add a text uicontrol to label the slider.
    txt = uicontrol('Style','text',...
        'Position',[400 45 120 20],...
        'String','Vertical Exaggeration');
    
    % Make figure visble after adding all components
    f.Visible = 'on';
    % This code uses dot notation to set properties. 
    % Dot notation runs in R2014b and later.
    % For R2014a and earlier: set(f,'Visible','on');
    

    function surfzlim(source,callbackdata)
        % initial display
        atlas_map=zeros(91,109,91);
        for ss=1:size(T,1) 
            atlas_map(F == ss) = T(ss,source.Value,1);
        end
        for ss = 1:size_range;
    
        nwe_B = B(:,:,range(ss));    
        new_F = atlas_map(:,:,range(ss));

        % To have a grayscale background, replicate image to 3-channels
        nwe_B = repmat(mat2gray(double(nwe_B),double(climB)),[1,1,3]);

        % Display the back image
        axes(ha(ss));
        hB = imagesc(nwe_B,climB);  axis image off;
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

end