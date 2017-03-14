function [ K_train_test_cn ] = kernel_preprocessing_test( K_train_test, K_test, train_process, ker_center )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

    % centering
    mean_K_train = train_process.mean;
    mean_K_train_test = mean(K_train_test, 1);
    
    K_train_test_c = center(K_train_test, mean_K_train, ker_center, mean_K_train', mean_K_train_test);
    K_test_c = center(K_test, mean_K_train, ker_center, mean_K_train_test', mean_K_train_test);
    
    % normalization
    K_train_test_cn = normmat(K_train_test_c, train_process.diag_c, diag(K_test_c));

end

