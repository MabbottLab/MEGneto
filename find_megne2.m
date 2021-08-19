function [ path ] = find_megne2(  )
%FIND_MEGNE2 Finds the path to the directory megne2 is running in
path = which(mfilename());
path = path(1:end-(length(['/' mfilename() '.m'])));
end

