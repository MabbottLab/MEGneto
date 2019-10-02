function [mymap_conn, mymap_diff] = create_colourmap(visual)

% create colourmap
cmap_list=brewermap('list');

tmp_conn = strsplit(visual.cmap_conn,'*'); % check if reversed colourmap
if length(tmp_conn) > 1
    cmap_conn =tmp_conn{1,2} ; 
else
    cmap_conn = visual.cmap_conn; 
end
if nnz(strcmp(cmap_list,cmap_conn)) > 0   % check if custom or standard colourmap
    mymap_conn = brewermap(100,visual.cmap_conn);
else
    mymap_conn = visual.cmap_conn;
end

tmp_diff = strsplit(visual.cmap_diff,'*'); % check if reversed colourmap
if length(tmp_diff) > 1
    cmap_diff = tmp_diff{1,2}; 
else
    cmap_diff = visual.cmap_diff; 
end
if nnz(strcmp(cmap_list,cmap_diff)) > 0 % check if custom or standard colourmap
    mymap_diff = brewermap(100,visual.cmap_diff);
else
    mymap_diff = visual.cmap_diff;
end


end