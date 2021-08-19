function [ conf ] = load_config( paths, file)
%LOAD_CONFIG Loads specified config file from the config directory, with
%some more verbose error handling specific to the pipeline
config_loc = [paths.conf_dir '/' file '.json'];
if ~exist(paths.anhome,'dir')
    error('Analysis path does not yet exist. Make sure to initialize with megne2setup')
end
if ~exist([paths.conf_dir '/' file '.json'],'file')
    error('Config json file for specified analysis step does not exist.')
end
conf = loadjson(config_loc);
conf = recursive_json_struct_string_to_func(conf);
end

