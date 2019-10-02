clear;
%%
% This script creates three files: In current directory (can move them
% later)
% 1) BinaryMatrix*.txt (ascii file)
% 2) BinaryMatrix*.edge (ascii file for BNV)
% 3) BinaryMatrix*.mat (matlab data file)
%   example output:
%           BinaryMatrix_cerebellum_patient_control_13-29Hz_2_25.mat
%           BinaryMatrix_<network>_<group1>_<group2>_<freq>_<threshold>.mat
%
% Input File : NBS saved file from GUI
%   example : NBS_cerebellum_patient_control_13-29Hz_2_25.mat 
%             NBS_<network>_<group1>_<group2>-<freq>_<threshold>.mat
% 
% 

%%
% User input : This needs to be changed

baseDir = '/data2/connectivity_study/nbs_analysis2/ventralAttention/NBS_results';



%%
%vector of frequences used from run_pilepline.m
%filt_freqs     = [2 3; 4 7; 8 12; 13 29; 30 59; 60 100; 101 150];


global nbs;

[LIST, ~] = glob([baseDir,'/NBS*.mat']);
nbsFiles = LIST;

for i=1%:length(nbsFiles)
    disp(nbsFiles{i})
    %load mat file 
    load(nbsFiles{i})
    
    temp_file = strsplit(nbsFiles{i},{'/'});
    temp_file2 = strsplit(temp_file{end},{'_','.'});
    
    adj = nbs.NBS.con_mat{1}+nbs.NBS.con_mat{1}';
    %save txt file
    filename_binary = ['BinaryMatrix_' temp_file2{2} '_' temp_file2{3} '_' temp_file2{4} '_' temp_file2{5} '_' temp_file2{6} '_' temp_file2{7} '_' temp_file2{8} '.txt'];
    fprintf('Saving txt file %s \n',filename_binary);
%     dlmwrite(filename_binary,full(adj),'delimiter',' ','precision','%d');
    %save edge file
    filename_edge = ['BinaryMatrix_' temp_file2{2} '_' temp_file2{3} '_' temp_file2{4} '_' temp_file2{5} '_' temp_file2{6} '_' temp_file2{7} '_' temp_file2{8} '.edge'];
    fprintf('Saving edge file %s \n',filename_edge);
%     dlmwrite(filename_edge,full(adj),'delimiter',' ','precision','%d');
    %save mat file
    filename_binarymat = ['BinaryMatrix_' temp_file2{2} '_' temp_file2{3} '_' temp_file2{4} '_' temp_file2{5} '_' temp_file2{6} '_' temp_file2{7} '_' temp_file2{8} '.mat'];
    fprintf('Saving mat file %s \n',filename_binarymat);
%     save(filename_binarymat,'adj')

end
