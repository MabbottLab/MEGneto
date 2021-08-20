function crossValidation(paths, num_interest, num_control)

config      = load_config(paths, paths.name);
config      = config.config;
band_names  = config.connectivity.freq_names;

for fq = 1:length(config.connectivity.filt_freqs)
    
    load([paths.anout_grp '/ml' band_names{fq} '_ml_data.mat']);
    
    % initialize variables
    for_later                         = struct();
    validationAccuracy                = [];
    finalPredictions                  = [];
    finalHyperparameters              = [];
    
    thresholds      = [];
    
    % Train-test 80/20 split
    total_num                         = num_interest + num_control;
    test_split                        = round(0.2 * total_num);
    split_1                           = round(test_split / 2);
    split_2                           = test_split - split_1;
    
    training_interest                 = num_interest - split_1;
    training_control                  = num_control - split_2;
    
    interest_test                     = randperm(num_interest, split_1);
    control_test                      = num_interest + randperm(num_control, split_2);
    training_data                     = ml_data;
    training_data(control_test, :)    = [];
    training_data(interest_test, :)   = [];
    testing_data                      = [ml_data(interest_test, :); ml_data(control_test, :)];
    
    % Initialize classes
    classes(1:num_interest)             = "interest";
    classes(num_interest + 1:total_num) = "control";
    training_classes(1:training_interest) = "interest";
    training_classes(training_interest + 1:training_interest + training_control) = "control";
    
    CV_SVM                            = struct();
    boxConstraint                     = [];
    kernelScale                       = [];
    list                              = [];
    temp                              = struct();
    
    % nested cross-validation to find optimal hyperparameters
    for ss = 1:height(training_data)
        innerCV         = training_data;
        innerCV(ss, :)  = [];
        CV_class        = training_classes;
        CV_class(ss)    = [];
        CV_test         = training_data(ss,:);
        CV_testclass    = classes(ss);
        % applied to (k-1) folds or groups from outer nest
        % optimize parameters, use to configure the model
        % evaluate on the remaining fold
        inner_c         = cvpartition(height(innerCV), 'KFold', height(innerCV));
        
        opts            = struct('Optimizer', 'bayesopt', 'ShowPlots', false, 'CVPartition', inner_c, ...
            'AcquisitionFunctionName', 'expected-improvement-plus');
        
        innerCV_SVM     = fitcsvm(innerCV, CV_class, 'KernelFunction', 'rbf', ...
            'OptimizeHyperparameters', 'auto', 'HyperparameterOptimizationOptions', opts, ...
            'ClassNames', {'interest'; 'control'});
        % validate on left-out observation
        
        temp_boxConstraint = innerCV_SVM.HyperparameterOptimizationResults.XAtMinObjective.BoxConstraint;
        temp_kernelScale   = innerCV_SVM.HyperparameterOptimizationResults.XAtMinObjective.KernelScale;
        
        [label,score] = predict(innerCV_SVM, CV_test);
        if label{1,1} == CV_testclass
            boxConstraint = [boxConstraint; temp_boxConstraint];
            kernelScale   = [kernelScale; temp_kernelScale];
        end
    end
    
    % thresholds for box constraint and kernel scale - can be changed
    % depending on sample size + other variables
    temp.hyperparameters = [boxConstraint kernelScale];
    tboxConstraint       = 15;
    tkernelScale         = 5;
    for ss = 1:length(temp.hyperparameters)
        if (temp.hyperparameters(ss, 1) < tboxConstraint) && (temp.hyperparameters(ss, 2) < tkernelScale)
            boxConstraint = temp.hyperparameters(ss, 1);
            kernelScale   = temp.hyperparameters(ss, 2);
            thresholds    = [thresholds; boxConstraint kernelScale];
        end
    end
    
    acc    = 0;
    for ss = 1:height(thresholds)
        classificationSVM_true = fitcsvm(...
            ml_data, classes, ...
            'KernelFunction', 'rbf', ...
            'PolynomialOrder', [], ...
            'KernelScale', thresholds(ss, 2), ...
            'BoxConstraint', thresholds(ss, 1), ...
            'Standardize', true, ...
            'ClassNames', {'interest'; 'control'});
        
        partitionedModel_true = crossval(classificationSVM_true, 'KFold', 20);
        
        % Compute validation predictions
        [validationPredictions_true, ~] = kfoldPredict(partitionedModel_true);
        
        % Compute validation accuracy
        validationAccuracy_true = 1 - kfoldLoss(partitionedModel_true, 'LossFun', 'ClassifError');
        
        % Test on reserved testing set
        [test_labels_true, test_accuracy_true] = predict(classificationSVM_true, testing_data);
        
        classificationSVM_false = fitcsvm(...
            ml_data, classes, ...
            'KernelFunction', 'rbf', ...
            'PolynomialOrder', [], ...
            'KernelScale', thresholds(ss, 2), ...
            'BoxConstraint', thresholds(ss, 1), ...
            'Standardize', false, ...
            'ClassNames', {'nonsocial'; 'social'});
        partitionedModel_false = crossval(classificationSVM_false, 'KFold', 20);
        
        % Compute validation predictions
        [validationPredictions_false, ~] = kfoldPredict(partitionedModel_false);
        
        % Compute validation accuracy
        validationAccuracy_false = 1 - kfoldLoss(partitionedModel_false, 'LossFun', 'ClassifError');
        
        % Test on reserved testing set
        [test_labels_false, test_accuracy_false] = predict(classificationSVM_false, testing_data);
        
        if (validationAccuracy_true > acc) || (validationAccuracy_false > acc)
            if validationAccuracy_true >= validationAccuracy_false
                %                     classificationSVM.(['perm' num2str(permutation)])   = classificationSVM_true;
                validationAccuracy                                  = validationAccuracy_true;
                acc                                                 = validationAccuracy_true;
                validationPredictions                               = validationPredictions_true;
                index                                               = ss;
                testingAccuracy                                     = test_accuracy_true;
                testingPredictions                                  = test_labels_true;
                standardized                                        = 'True';
            else
                %                     classificationSVM.(['perm' num2str(permutation)])   = classificationSVM_false;
                validationAccuracy                                  = validationAccuracy_false;
                acc                                                 = validationAccuracy_false;
                validationPredictions                               = validationPredictions_false;
                index                                               = ss;
                testingAccuracy                                     = test_accuracy_false;
                testingPredictions                                  = test_labels_false;
                standardized                                        = 'False';
            end
            %             else
            %                 validationAccuracy                             = [validationAccuracy; 0];
            %                 finalPredictions                               = [finalPredictions; num2cell(0)];
        end
    end
    for_later.validationAccuracy     = validationAccuracy;
    for_later.validationLabels       = validationPredictions;
    for_later.usedthreshold          = thresholds(index, :);
    for_later.thresholds             = thresholds;
    for_later.testingAccuracy        = testingAccuracy;
    for_later.testingLabels          = testingPredictions;
    for_later.standardized           = standardized;
    save([paths.anout_ML '/' band_names{fq} 'MLresults.mat'], 'for_later','-v7.3');
end
end

