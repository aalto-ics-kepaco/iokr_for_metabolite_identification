function [ KX_train, KX_train_test ] = mkl_combine_train_test( KX_list, train_set, test_set, w, ker_center)
%======================================================
% DESCRIPTION:
% Kernel combination using the weights learnt with some MKL algorithm
%
% INPUTS:
% KX_list:      cell array containing the Gram matrices of the different input kernels (size n_kx*1)
% train_set:    vector containing the training indices (length n_train)
% test_set:     vector containing the test indices (length n_test)
% w:            vector containing the weights learned with the 
%               mkl_weight.m function (length: n_kx)
% ker_center:   binary value indicating if the base kernels should be centered
%
% OUTPUTS:
% KX_train:     combined kernel matrix on the training set (size n_train*n_train)
% KX_train_test:combined kernel matrix between the training set and the
%               test set (size n_train*n_test)
%
%======================================================

    n_kx = length(KX_list); % number of input kernels
    n_train = length(train_set);
    n_test = length(test_set);

    % Combination of the kernels
    KX_train = zeros(n_train,n_train);
    KX_train_test = zeros(n_train,n_test);
    for i = 1:n_kx
        % Centering and normalization
        [KX_train_nc, KX_train_test_nc] = input_kernel_center_norm( KX_list{i}, train_set, test_set, ker_center);
        KX_train = KX_train + w(i) * KX_train_nc;
        KX_train_test = KX_train_test + w(i) * KX_train_test_nc;
    end

end

