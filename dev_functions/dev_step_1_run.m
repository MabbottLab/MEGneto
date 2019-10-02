% copyfile('/mnt/sda/ming/Summer2019/code/MEGne2/configs/config.json', ...
%     '/mnt/sda/ming/Summer2019/code/MEGne2/templates/test_1.json')
analyses = ["EL","TP"];
project_path = '/mnt/sda/ming/Summer2019';
rawdata_path = '/mnt/sda/ming/Summer2019/rawdata';
mri_path = '/mnt/sda/ming/Summer2019/rawdata/MRIs';
for ii = analyses
    jj = char(ii);
    analysis_name = jj;
    paths = megne2setup(project_path, analysis_name, rawdata_path, mri_path, false);
    fcp_2_PreprocessingICA(paths, false, true)
end