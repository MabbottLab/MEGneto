function ml_datamaker(paths, regions_of_interest, test, alpha, direction)

%% Initialize
config      = load_config(paths, paths.name);
config      = config.config;
band_names  = config.connectivity.freq_names;
ROIs        = config.connectivity.ROIs;
flag        = 0;
reverse     = 0;
indices     = regions_of_interest;

if isempty(regions_of_interest)
    ROIs = config.connectivity.ROIs;
else
    ROIs = regions_of_interest;
end

if strcmp(test, 'one-tailed')
    threshold = alpha;
else
    threshold = alpha/2;
    flag = 1;
end

if strcmp(direction, 'smaller')
    threshold = 1-alpha;
    reverse = 1;
end
   

for fq = 1:length(config.connectivity.filt_freqs)
    load([paths.anout_grp '/fcp_5_interest_allParticipants_conn_mats_' config.connectivity.method '.mat']);
    interest = all_conn_mat(ROIs,ROIs,:,fq);
    load([paths.anout_grp '/fcp_5_control_allParticipants_conn_mats_' config.connectivity.method '.mat']);
    control = all_conn_mat(ROIs,ROIs,:,fq);
    
    load([paths.anout_grp '/permstats' band_names{fq} '_p_val.mat']);
    indices = [];
    if reverse == 1 && flag == 1
        fprintf("Not possible");
    elseif flag == 1
        indices = find(p_val < threshold);
        indices = [indices; find(p_val > 1-threshold)];
    else
        indices = find(p_val < threshold);
    end
    temp = [];
    temp(:, 2)  = mod(indices, length(ROIs));
    for i = 1:height(temp)
        if temp(i,2) > 0
            temp(i,1) = floor(indices(i)/length(ROIs)) + 1;
        else
            temp(i,1) = floor(indices(i)/length(ROIs));
        end
        if temp(i,2) == 0
            temp(i,2) = length(ROIs);
        end
    end
    
    for i = 1:height(temp)
        interest_ml(:,i) = squeeze(interest(temp(i,1), temp(i,2), :));
        control_ml(:, i) = squeeze(control(temp(i,1), temp(i,2), :));
    end
    ml_data = [];
    temp = [];
    ml_data = [interest_ml; control_ml];
    temp = unique(ml_data', 'rows');
    ml_data = temp';
    save([paths.anout_grp '/ml' band_names{fq} '_ml_data.mat'], 'ml_data','-v7.3');
end