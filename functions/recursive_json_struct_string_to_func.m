function y = recursive_json_struct_string_to_func(y)
z = fieldnames(y);
functerprint = {'function','type','file','workspace','within_file_path'};
if all(ismember(functerprint,z)) && length(z) == length(functerprint)
    y = str2func(y.function);
    return
end
for ii = 1:length(z)
    switch class(y.(z{ii}))
        case 'struct'
            y.(z{ii}) = recursive_json_struct_string_to_func( y.(z{ii}) );
    end
end
end