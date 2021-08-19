function [centroid] = find_centroid(aal_node,sourcemodel)
% aal_node: array comprised of indices of sources in region
%       example: aal_node=find(source_aal.tissue==1);
% sourcemodel: is the source model structure

if length(aal_node)==1 % only one source
    centroid=aal_node; 
else % find centroid by computing minimum average distance between sources
    dist_array=zeros(size(aal_node,1),1);
    stacked_avg_dist=zeros(size(aal_node,1),1);
    for j=1:size(aal_node,1) % each source becomes the reference node
        ref_node=aal_node(j); 
        ref_coords=sourcemodel.pos(ref_node,:);
        for k=1:size(aal_node,1) % determine distance between reference node and every other source
            diff_vector=sourcemodel.pos(aal_node(k),:)-ref_coords(1,:);
            dist_array(k,1)=(diff_vector(1))^2+(diff_vector(2))^2+(diff_vector(3))^2; % magnitude squared
        end
        avg_dist=mean(dist_array,1); % compute average distance for reference node
        stacked_avg_dist(j,1)=avg_dist; % stack average distance for each reference node
    end
    [~,I]=min(stacked_avg_dist); % find reference node with minimum average distance
    centroid=aal_node(I); % that becomes the centroid
end


return