function [ status ] = encodeH264( target, varargin )
%encodeH264 - Transcode the .avi video given in `target` into .mp4
%             
% Requires ffmpeg.
% 
% Simeon Wong
% 2014 June 18

%% Parse Parameters
p = inputParser;

addOptional(p, 'qp', 25);
addOptional(p, 'savepath', []);

parse(p, varargin{:});

if isempty(p.Results.savepath)
    [path,name,~] = fileparts(target);
    save_path = fullfile(path, [name, '.mp4']);
else
    save_path = p.Results.savepath;
end

%% Encode

[status,output] = system(sprintf('ffmpeg -i %s -vcodec libx264 -qp %d -y %s', target, p.Results.qp, save_path));
if status == 0
    delete(target);
else
    disp(output);
    warning('Error encoding video!');
end


end

