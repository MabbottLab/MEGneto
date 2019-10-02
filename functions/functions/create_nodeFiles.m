clear;
%%
% This script creates ten files: In current directory (can move them
% later)
% 1) significat network for group1 (matlab data file)
%          example: patients_group_13-29Hz_signetwork_cerebellum.mat
%                   patients_group_<freq>_signetwork_<network>.mat
% 2) significat network for group2 (matlab data file) 
%          example: controls_group_13-29Hz_signetwork_cerebellum.mat
%                   controls_group_<freq>_signetwork_<network>.mat
% 3) strength files for group1 (matlab data file)
%          example: patients_group_13-29Hz_cerebellum_strength.mat
%                   patients_group_<freq>_<network>_stength.mat
% 4) strength files for group2 (matlab data file)
%          example: controls_group_13-29Hz_cerebellum_strength.mat
%                   controls_group_<freq>_<network>_stength.mat
% 5) node file: for group1 (ascii file for BNV) with spaces separating
% columns (if you need tabs you will need to change the code)
%          example: cerebellum_patient_group_13-29Hz_nodes.txt
%                   <network>_patient_group_<freq>_nodes.txt
% 6) node file: for group1 (ascii file for BNV)
%          example: cerebellum_patient_group_13-29Hz_nodes.node
%                   <network>_patient_group_<freq>_nodes.node
% 7) node file: for group2 (ascii file for BNV)
%          example: cerebellum_control_group_13-29Hz_nodes.txt
%                   <network>_control_group_<freq>_nodes.txt
% 8) node file: for group2 (ascii file for BNV)
%          example: cerebellum_control_group_13-29Hz_nodes.node
%                   <network>_control_group_<freq>_nodes.node
% 9) node file: for difference between group 1 & group2 (ascii file for BNV)
%          example: cerebellum_control_group_13-29Hz_nodes.txt
%                   <network>_control_group_<freq>_nodes.txt
% 10) node file: for difference between group 1 & group2 (ascii file for BNV)
%          example: cerebellum_difference_control_group_13-29Hz_nodes.node
%                   <network>_difference_control_group_<freq>_nodes.node
%
%
% Input Files : 
% 1) NBS saved file from fcp_5_NBS.mat (stored in network directory)
%   example : NBS_13-29Hz.mat 
%             NBS_<freq>.mat
% 2) BinaryMatrix 
%   example: BinaryMatrix_cerebellum_patient_control_13_29Hz_2_25.mat
%            BinaryMatrix_<network>_<group1>_<group2>_<freq>_<threshold>.mat
%
% 3) Template Node files for 6 networks (saved in main MEGneto directory)
%   example: cerebellum_template_node_file.txt
%
%%

% addpath for BCT toolbox
addpath('/data2/connectivity_study/analysis_tools/BCT/2015_01_25_BCT');
% addpath('/home/sonya/matlab/BCT');

%%
% User input : This needs to be changed
% Setup main variables
%Where NBS_8-12Hz.mat for that network is stored
baseDir = '/data2/connectivity_study/nbs_analysis2/visualNetwork';
%Where binary matrix files are stored
binaryDir = '/data2/connectivity_study/nbs_analysis2/visualNetwork/binary_matrix';
%Where the template nodefiles are
templateNodeDir = '/data2/connectivity_study/analysis_tools/MEGneto/functions/NodeFiles';

%number of participants in each group
numPs  = 24;
%current network : in baseDir and binaryDir
network = 'visual';
%vector of frequences used from run_pilepline.m
% filt_freqs     = [2 3; 4 7; 8 12; 13 29; 30 59; 60 100; 101 150]; %do not change this
filt_freqs     = [101 150; 13 29; 2 3; 30 59; 4 7; 60 100; 8 12];
%significant frequencies of the network
sig_filt_freq = [1,0,0,0,0,0,0]; %0 = non-significant; 1 = significant

%%
% Path to subjects directories
[LIST, ~] = glob([baseDir,'/NBS*Hz.mat']);
nbsFiles = LIST;

[LIST, ~] = glob([binaryDir,'/BinaryMatrix*.mat']);
binaryFiles = LIST;

[LIST, ISDIR] = glob([templateNodeDir,'/*.txt']);
templateNodeFiles = LIST;

for i=1:length(nbsFiles) %Fix!!!!
    load(nbsFiles{i})
    disp(nbsFiles{i})
    temp_freq = strsplit(nbsFiles{i},{'/','_','.'});
    temp_freq2 = temp_freq{end-1};
    ch_freq = strsplit(temp_freq2,'-');
    ind = find(filt_freqs == str2num(ch_freq{1}));
    if sig_filt_freq(ind) == 1
        fprintf('Workign on frequency : %s \n',temp_freq2)
        %loads file called data_matrix contails both groups
        %seperate between the two groups
        group1 = data_matrix(:,:,1:numPs);
        group2 = data_matrix(:,:,numPs+1:end);
        mean_group1 = mean(group1,3);
        mean_group2 = mean(group2,3);

        %Load binary matrix file
        %ind2 = strfind(binaryFiles,'beta'); %for testing FIX!!!!
        ind2 = strfind(binaryFiles,temp_freq2);
        ind3 = find(not(cellfun('isempty',ind2)));
        fprintf('Loading binary file : %s \n',binaryFiles{ind3});
        load(binaryFiles{ind3}); %loads variable adj
        sig_network_group1 = times(mean_group1,adj);
        sig_network_group2 = times(mean_group2,adj);
        %
        fileName_group1 = ['patients_group_' temp_freq2 '_signetwork_' network '.mat'];
        fileName_group2 = ['controls_group_' temp_freq2 '_signetwork_' network '.mat'];
        save(fileName_group1,'sig_network_group1')
        save(fileName_group2,'sig_network_group2')
        
        %Now lets calculate strenghts :)
        fprintf('Calculating strenghts for group: %s \n', 'patients');
        strength_group1 = strengths_und(sig_network_group1);      
        fileName_group1 = ['patients_group_' temp_freq2 '_' network '_strength.mat'];
        save(fileName_group1,'strength_group1')
        
        fprintf('Calculating strenghts for group: %s \n', 'controls');
        strength_group2 = strengths_und(sig_network_group2);
        fileName_group2 = ['controls_group_' temp_freq2 '_' network '_strength.mat'];
        save(fileName_group2,'strength_group2')
        
        %now for the node files
        %1) load/read in base file to modify
        ind2 = strfind(templateNodeFiles,network); %
        ind3 = find(not(cellfun('isempty',ind2)));
        fprintf('Loading template node file : %s \n',templateNodeFiles{ind3});
        fileID = fopen(templateNodeFiles{ind3},'r');
        A = textscan(fileID, '%d %d %d %d %f %s'); %
        fclose(fileID);
        Abkp = A; %backup the default node file
        
        %for group1 : node file
        temp_strength_group1 = full(strength_group1); %
        A{1,5} = A{1,5}+temp_strength_group1';
        fileName = [network '_patient_group_' temp_freq2 '_nodes.txt' ];
        fileID = fopen(fileName,'w');
        for ii=1:length(A{1,1})
            %fprintf(fileID, '%d\t%d\t%d\t%d\t%f\t%s \n', A{1,1}(i),
            %A{1,2}(i), A{1,3}(i), A{1,4}(i), A{1,5}(i), A{1,6}{i}); %with
            %tab
            fprintf(fileID, '%d %d %d %d %f %s \n', A{1,1}(ii), A{1,2}(ii), A{1,3}(ii), A{1,4}(ii), A{1,5}(ii), A{1,6}{ii});
        end    
        fclose(fileID);
        temp_fileName = strsplit(fileName,{'.'});
        command = ['cp ' fileName ' ' temp_fileName{1} '.node']; %fileName ' 'fileName];
        [status,cmdout] = system(command);
        A = Abkp;
        %%
        %for group2
        temp_strength_group2 = full(strength_group2); %test
        A{1,5} = A{1,5}+temp_strength_group2';
        fileName = [network '_control_group_' temp_freq2 '_nodes.txt' ];
        fileID = fopen(fileName,'w');        
        for ii=1:length(A{1,1})
            %fprintf(fileID, '%d\t%d\t%d\t%d\t%f\t%s \n', A{1,1}(i),
            %A{1,2}(i), A{1,3}(i), A{1,4}(i), A{1,5}(i), A{1,6}{i}); %with
            %tab
            fprintf(fileID, '%d %d %d %d %f %s \n', A{1,1}(ii), A{1,2}(ii), A{1,3}(ii), A{1,4}(ii), A{1,5}(ii), A{1,6}{ii});
        end    
        fclose(fileID);
        temp_fileName = strsplit(fileName,{'.'});
        command = ['cp ' fileName ' ' temp_fileName{1} '.node']; %fileName ' 'fileName];
        [status,cmdout] = system(command);
        A = Abkp;
        %for group1-group2
        temp_strength_diff = temp_strength_group1-temp_strength_group2; %test
        A{1,5} = A{1,5}+temp_strength_diff';
        fileName = [network '_difference_patient_control_group_' temp_freq2 '_nodes.txt' ];
        fileID = fopen(fileName,'w');        
        for ii=1:length(A{1,1})
            %fprintf(fileID, '%d\t%d\t%d\t%d\t%f\t%s \n', A{1,1}(i),
            %A{1,2}(i), A{1,3}(i), A{1,4}(i), A{1,5}(i), A{1,6}{i}); %with
            %tab
            fprintf(fileID, '%d %d %d %d %f %s \n', A{1,1}(ii), A{1,2}(ii), A{1,3}(ii), A{1,4}(ii), A{1,5}(ii), A{1,6}{ii});
        end    
        fclose(fileID);
        temp_fileName = strsplit(fileName,{'.'});
        command = ['cp ' fileName ' ' temp_fileName{1} '.node']; %fileName ' 'fileName];
        [status,cmdout] = system(command);
    else 
         fprintf('Frequency : %s is not significant \n',temp_freq2)
    end
      
end