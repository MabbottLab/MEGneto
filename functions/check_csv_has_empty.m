function out = check_csv_has_empty(input_file)
%CHECK_CSV_HAS_EMPTY accepts an input path leading to a
%readtable-compatible file and checks if it contains any empty values
try
    tab = table2cell(readtable(input_file));
    if isempty(tab)
        out = 1;
    else
        out = sum(any(cellfun('isempty',tab)));
    end
catch
    out = 1;
end
