function [normalized_weigthing] = region_source_weighting(centroid,aal_node,sourcemodel)
% aal_node: array comprised of indices of sources in region
%       example: aal_node=find(source_aal.tissue==1);
% sourcemodel: is the source model structure
% centroid: source index that represents the centre source

centroid_dist=zeros(size(aal_node,1),1);
centroid_coords=sourcemodel.pos(centroid,:);
for k=1:size(aal_node,1) % determine distance between centroid and every other source
    cent_diff_vector=sourcemodel.pos(aal_node(k),:)-centroid_coords(1,:);
    centroid_dist(k,1)=(cent_diff_vector(1))^2+(cent_diff_vector(2))^2+(cent_diff_vector(3))^2; % magnitude squared
end
%                 a=1; % amplitude
%                 b=0; % mean
%                 c=1; % variance in cm
%                 x=centroid_dist; % distance b/w centroid and other sourced in cm
%                 weighting=a*exp(-(x-b).^2/(2*c^2)); % gaussian function
max_dist=max(centroid_dist); % in cm
reversed_dist=max_dist-centroid_dist;
weighting=reversed_dist.^2; % increase contribution from central voxels
normalized_weigthing=weighting/sum(weighting); % normalize weighting to equal to 1

normalized_weigthing(isnan(normalized_weigthing)) = 1; % if only one source in region give it a weighting of one
 
return