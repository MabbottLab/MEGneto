function [labels,output_pls,sig_nodes] = plot_matrix_connectivity_AAL(fc,thr,thrpar,figTitle,regions)


load('AAL116_region_names'); % 'region_labels'
nchan = length(sort(cell2mat(regions)));
regions=sort(cell2mat(regions));


% make matrix as output only (not visualization)----------------------------------------------------

%region names
roi_names=region_labels(regions)';

% thresholding to get network
CC = fc;
if strcmp(thrpar,'upper') == 1,
    CC(CC<thr) = 0;
    %CC(CC>=thr) = 1;
    %my_color = 'r';
elseif strcmp(thrpar,'lower') == 1,
    CC(CC>thr) = 0;
    %CC(CC<=thr) = -1;
    %my_color = 'b';
elseif strcmp(thrpar,'all') == 1,
    CC(CC>-thr & CC<thr)=0;
    disp('doing nothing');
else
    disp('Wrong parameter: ''lower'' or ''upper'' ');
    return
end   

% function outputs
output_pls=CC;
sig_nodes=output_pls;
sig_nodes(find(sig_nodes~=0))=1;
labels=roi_names;

% make matrix for visualization--------------------------------------------

% separates left and right regions
right=find(mod(regions,2)==0);
left=find(mod(regions,2)~=0);
i_left_right = [regions(left) regions(fliplr(right))]; 
% region names
roi_names = region_labels(i_left_right)';

% rearange fc into left and right
fc_right=find(mod(1:nchan,2)==0);
fc_left=find(mod(1:nchan,2)~=0);
fc_order=[fc_left fliplr(fc_right)]; 

trans_mat = zeros(nchan); 
trans_mat(:,1:nchan) = fc(:,fc_order); 
CC=zeros(nchan);
CC(1:nchan,:)=trans_mat(fc_order,:);

% thresholding to get network
if strcmp(thrpar,'upper') == 1,
    CC(CC<thr) = 0;
    %CC(CC>=thr) = 1;
    %my_color = 'r';
elseif strcmp(thrpar,'lower') == 1,
    CC(CC>thr) = 0;
    %CC(CC<=thr) = -1;
    %my_color = 'b';
elseif strcmp(thrpar,'all') == 1,
    CC(CC>-thr & CC<thr)=0;
    disp('doing nothing');
else
    disp('Wrong parameter: ''lower'' or ''upper'' ');
    return
end

% VISUALIZATION------------------------------------------------------------

% make nice labels for figure
clear a b;
for jj=1:nchan, % shorten region names
  [a,b] = strtok(roi_names{jj},'_');
  Labels{jj} = [upper(a)];
end    

adjMat = CC;

% figure title
figure; 
set(gcf,'color',[1 1 1],'Name',figTitle,'renderer','zbuffer');
set(gcf, 'Units', 'centimeters','PaperPositionMode','auto');
set(gcf, 'Position', [1 1 31 31]);

subplot(100,100,[8 10000]);
adjMat = [[adjMat,zeros(size(adjMat,1),1)]; zeros(1,size(adjMat,1)+1)]; %pads matrix to be 117x117
h = pcolor(adjMat);

% define colormap
my_map = bipolar;
my_map = colormap;
if strcmp(thrpar,'upper') == 1,
        my_map = my_map(33:64,:);
        % my_map = colormap;
        c2 = max(max(abs(adjMat)));
        c1 = 0;
elseif strcmp(thrpar,'lower') == 1,
        my_map = my_map(1:32,:);
        c1 = min(min(adjMat));
        c2 = 0;
elseif strcmp(thrpar,'all') == 1, 
    c2 = max(max(abs(adjMat)));
    c1 = -c2;
    my_map = colormap;
else
    disp('Wrong parameter: ''lower'' or ''upper'' ');
    return
end
colormap (my_map);
caxis([c1 c2]);

% x and y axis labels
set(gca,'TickLength',[0 0],'TickDir','in');
set(gca, 'FontSize', 6,'FontWeight','bold');
set(gca,'YTick',[1:numel(Labels)]+0.6,'YTickLabel',Labels);
set(gca,'XTick',[1:numel(Labels)]+0.6,'XTickLabel',Labels);
xticklabel_rotate([],90,[],'FontSize', 6,'FontWeight','bold');
pos=get(gca,'position');

% % vertical lines
% annotation('line',[pos(1)+0.96*0.502*pos(3) pos(1)+0.96*0.502*pos(3)],[pos(2) pos(4)+pos(2)],'linestyle','-','color','k','linewidth',2.2);
% annotation('line',[pos(1)+0.96*16.02/90*pos(3) pos(1)+0.96*16.02/90*pos(3)],[pos(2) pos(4)+pos(2)],'linestyle','-','color','g','linewidth',1.8);
% annotation('line',[pos(1)+0.96*24.02/90*pos(3) pos(1)+0.96*24.02/90*pos(3)],[pos(2) pos(4)+pos(2)],'linestyle','-','color','g','linewidth',1.8);
% annotation('line',[pos(1)+0.96*30.01/90*pos(3) pos(1)+0.96*30.01/90*pos(3)],[pos(2) pos(4)+pos(2)],'linestyle','-','color','g','linewidth',1.8);
% annotation('line',[pos(1)+0.96*39.02/90*pos(3) pos(1)+0.96*39.02/90*pos(3)],[pos(2) pos(4)+pos(2)],'linestyle','-','color','g','linewidth',1.8);
% annotation('line',[pos(1)+0.96*(16.02+45)/90*pos(3) pos(1)+0.96*(16.02+45)/90*pos(3)],[pos(2) pos(4)+pos(2)],'linestyle','-','color','g','linewidth',1.8);
% annotation('line',[pos(1)+0.96*(24.02+45)/90*pos(3) pos(1)+0.96*(24.02+45)/90*pos(3)],[pos(2) pos(4)+pos(2)],'linestyle','-','color','g','linewidth',1.8);
% annotation('line',[pos(1)+0.96*(30.01+45)/90*pos(3) pos(1)+0.96*(30.01+45)/90*pos(3)],[pos(2) pos(4)+pos(2)],'linestyle','-','color','g','linewidth',1.8);
% annotation('line',[pos(1)+0.96*(39.02+45)/90*pos(3) pos(1)+0.96*(39.02+45)/90*pos(3)],[pos(2) pos(4)+pos(2)],'linestyle','-','color','g','linewidth',1.8);
% % horizontal lines
% annotation('line',[pos(1) pos(1)+pos(3)*0.96],[pos(2)+0.50*pos(4) pos(2)+0.50*pos(4)],'linestyle','-','color','k','linewidth',2.2);
% annotation('line',[pos(1) pos(1)+pos(3)*0.96],[pos(2)+16/90*pos(4) pos(2)+16/90*pos(4)],'linestyle','-','color','g','linewidth',1.8);
% annotation('line',[pos(1) pos(1)+pos(3)*0.96],[pos(2)+24/90*pos(4) pos(2)+24/90*pos(4)],'linestyle','-','color','g','linewidth',1.8);
% annotation('line',[pos(1) pos(1)+pos(3)*0.96],[pos(2)+30/90*pos(4) pos(2)+30/90*pos(4)],'linestyle','-','color','g','linewidth',1.8);
% annotation('line',[pos(1) pos(1)+pos(3)*0.96],[pos(2)+39/90*pos(4) pos(2)+39/90*pos(4)],'linestyle','-','color','g','linewidth',1.8);
% annotation('line',[pos(1) pos(1)+pos(3)*0.96],[pos(2)+(16+45)/90*pos(4) pos(2)+(16+45)/90*pos(4)],'linestyle','-','color','g','linewidth',1.8);
% annotation('line',[pos(1) pos(1)+pos(3)*0.96],[pos(2)+(24+45)/90*pos(4) pos(2)+(24+45)/90*pos(4)],'linestyle','-','color','g','linewidth',1.8);
% annotation('line',[pos(1) pos(1)+pos(3)*0.96],[pos(2)+(30+45)/90*pos(4) pos(2)+(30+45)/90*pos(4)],'linestyle','-','color','g','linewidth',1.8);
% annotation('line',[pos(1) pos(1)+pos(3)*0.96],[pos(2)+(39+45)/90*pos(4) pos(2)+(39+45)/90*pos(4)],'linestyle','-','color','g','linewidth',1.8);

% % vertical lines
annotation('line',[pos(1)+0.96*58.02/116*pos(3) pos(1)+0.96*58.02/116*pos(3)],[pos(2) pos(4)+pos(2)],'linestyle','-','color','k','linewidth',2.2);
% % horizontal lines
annotation('line',[pos(1) pos(1)+pos(3)*0.96],[pos(2)+58/116*pos(4) pos(2)+58/116*pos(4)],'linestyle','-','color','k','linewidth',2.2);

% x and y titles
title(figTitle,'FontSize',30,'FontWeight','bold');
xlabel('LEFT                                          RIGHT','FontSize',20,'FontWeight','bold');
ylabel('LEFT                                          RIGHT','FontSize',20,'FontWeight','bold');
set(gca,'pos',[pos(1) pos(2) pos(3)*0.96 pos(4)]);
h = colorbar('location','EastOutside','position',[pos(1)+pos(3) pos(2)+pos(4)*0.4 0.02 pos(4)*0.3]);
set(h,'Fontsize',12,'FontWeight','bold');

return

