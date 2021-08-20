function permutation_statistics(paths, regions_of_interest)

config      = load_config(paths, paths.name);
config      = config.config;
band_names  = config.connectivity.freq_names;


for fq = 1:length(config.connectivity.filt_freqs)
    fprintf('Analyzing the %s band!\n', band_names{fq});
    if isempty(regions_of_interest)
        indices = 1:length(config.connectivity.ROIs);
    else
        indices = regions_of_interest;
    end
    
    num_unique  = (length(indices) * (length(indices) - 1) / 2);    % number of unique connections
    
    load([paths.anout_grp '/fcp_5_' band_names{fq} '_permutations_mats_' config.connectivity.method '.mat']);
    wPLI_control = all_conn_mat(indices, indices, :, :);
    
    wPLI_control = wPLI_control(:,:,:,1:10);
    
    load([paths.anout_grp '/fcp_5_control_allParticipants_conn_mats_' config.connectivity.method '.mat']);
    control_ttest = squeeze(all_conn_mat(indices, indices, :, fq));
    
    control_ttest = control_ttest(:,:,1:10);
    
    load([paths.anout_grp '/fcp_5_interest_allParticipants_conn_mats_' config.connectivity.method '.mat']);
    wPLI_interest = all_conn_mat(indices, indices, :, fq);
    
    wPLI_interest = wPLI_interest(:,:,1:10);
    
    social_t    = [];
    control_t   = [];
    t_stat      = [];
    tmp         = [];
    comparisons = [];
    temp        = [];
    C           = [];
    tmppval     = [];
    
    n = 1;
    for i = 1:length(indices) - 1
        for j = i+1:length(indices)
            interest_t(:, n) = squeeze(wPLI_interest(i, j, :));
            n = n + 1;
        end
    end
    
    n = 1;
    for i = 1:length(indices) - 1
        for j = i+1:length(indices)
            control_t(:, n) = squeeze(control_ttest(i, j, :));
            n = n + 1;
        end
    end
    
    for i = 1:num_unique
        [~, ~, ~, stats] = ttest2(interest_t(:, i), control_t(:,i));
        tmp(1, i) = stats.tstat;
    end
    t_stat = tmp;
    
    
    for jj = 1:500
        n = 1;
        for i = 1:length(indices) - 1
            for j = i+1:length(indices)
                comparisons(:,n,jj) = squeeze(wPLI_control(i,j,jj,:));
                n = n + 1;
            end
        end
    end
    
    for ii = 1:6000
        temp(ii, 1:2) = randperm(500,2);
    end
    C = unique(temp, 'rows');
    temp = C;
    
    for ii = 1:length(temp)
        tmp1 = comparisons(:,:,temp(ii,1));
        tmp2 = comparisons(:,:,temp(ii,2));
        for jj = 1:num_unique
            [~, ~, ~, stats] = ttest2(tmp1(:,jj), tmp2(:,jj));
            control_stat(jj, ii) = stats.tstat;
        end
    end
    
    % p-val for each connection
    for nn = 1:num_unique
        tmp = control_stat(nn, :);
        tmp = tmp(tmp > t_stat(nn));
        tmppval(nn) = length(tmp)/length(temp);
    end
    
    % mean difference in wPLI between conditions
    for i = 1:num_unique
        diff(i) = mean(interest_t(:,i) - control_t(:,i), 1);
    end
    
    n = 1;
    
    for i = 1:length(indices) - 1
        p_val(i, i) = NaN;
        p_val(length(indices), length(indices)) = NaN;
        for j = i+1:length(indices)
            p_val(j,i) = tmppval(n);
            p_val(i,j) = tmppval(n);
            n = n+1;
        end
    end
    save([paths.anout_grp '/permstats' band_names{fq} '_p_val.mat'], 'p_val','-v7.3');
    save([paths.anout_grp '/permstats' band_names{fq} '_results.mat'], 'tmppval', 't_stat', 'diff', '-v7.3'); 
    
end
right_now = clock;
fprintf('%02.f:%02.f:%02.f ============== Finished Processing ====================\n', ...
    right_now(4:6))
end
    
    