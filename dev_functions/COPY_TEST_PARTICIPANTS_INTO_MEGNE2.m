%COPY_TEST_PARTICIPANTS
csvs = glob([paths.conf_dir '/subj_fcp*.csv']);
for ii = 1:length(csvs)
    evalstring = ['cp /mnt/sda/ming/Summer2019/notes/' lower(paths.name) '_dsfiles.csv ' csvs{ii}];
    system(evalstring);
end