%% Getting individual connectivity values

subj1 = {'POGO_11' 'POGO_13' 'POGO_23' 'POGO_25' 'POGO_29' 'POGO_31' 'POGO_32' 'POGO_34' 'ST02' 'ST03' ...
      'ST05' 'ST06' 'ST08' 'ST14' 'ST20' 'ST27' 'ST36' 'ST45' ...
      'ST47' 'ST48' 'ST49' 'ST50' 'ST51' 'ST56'};
  
subj2 = {'H_02_48' 'H_03_49' 'H_05_50' 'ST10' 'ST12' 'ST17' 'ST18' 'ST30' 'ST33' 'ST40' ...
    'T_03_19' 'T_04_14' 'T_08_35' 'T_09_37' 'T_11_40' 'T_13_42' 'T_14_S05' 'T_15_43' 'T_17_44' 'T_19_46' ...
    'T_20_52' 'T_21_53' 'T_24_55' 'T_27_59'};

subj = [subj2,subj1];

% projdir = '/data2/connectivity_study/ANALYSIS/controls_patients_pli/';
projdir = '/data2/connectivity_study/nbs_analysis2/languageNetwork/';
file_in = strcat(projdir, '/NBS_13-29Hz.mat');
% 
% file_in = @(id) strcat(projdir, '/', id, '/ANALYSIS/conn_', id, '_116pli.mat');

inFile = load(file_in, 'data_matrix'); %matrix (116x116x60)
% maskB = load('/data2/connectivity_study/nbs_analysis2/visualNetwork/binary_matrix/BinaryMatrix_visual_PatientControl_101-150Hz_thr_2_23.mat'); 
tempFile = struct2array(inFile);
all_sbj = zeros(6,48);
all_temp_mask_beta = zeros(6,6,48);

for ss = 1:length(subj)
   
    temp_beta = tempFile(:,:,ss);
%     temp_mask_beta = temp_beta.*maskB.adj;
    all_temp_mask_beta(:,:,ss) = temp_beta;%temp_mask_beta;
    [str] = strengths_und(temp_beta);

%     str = sum(temp_mask_beta);
    disp(str);
    all_sbj(:,ss) = str;

end

% save('ind_node_str_language_4-7Hz.mat','all_sbj');
% save('ind_connectivity_language_4-7Hz_mask.mat','all_temp_mask_beta');

