function [ K_train_nc, K_train_test_nc] = input_kernel_center_norm( K, train_set, test_set, ker_center)
%======================================================
% DESCRIPTION:
% Centering and normalization of a Gram matrix

% INPUTS
% K:            Gram matrix
% train_set:    vector containing the indices of the training examples
% test_set:     vector containing the indices of the test_examples
% ker_center:   binary value indicating if the Gram matrix has to be centered or not
%
% OUTPUTS:
% K_train_nc:   Gram matrix between training examples
% K_train_test_nc: Gram matrix between training and test examples
%
%======================================================

    K_train = K(train_set,train_set);
    K_test = K(test_set,test_set);
    K_train_test = K(train_set,test_set);
    
    mean_K_train = mean(K_train,1);
    mean_K_train_test = mean(K_train_test,1);

    if ker_center == 1
        % Centering
        K_train_c = center(K_train, mean_K_train, ker_center);
        K_test_c = center(K_test, mean_K_train, ker_center, mean_K_train_test', mean_K_train_test);
        K_train_test_c = center(K_train_test, mean_K_train, ker_center, mean_K_train', mean_K_train_test);

        % Normalization
        K_train_nc = normmat(K_train_c);
        K_train_test_nc = K_train_test_c ./ sqrt(diag(K_train_c) * diag(K_test_c)');
        
    else % because input kernels are already normalized in our application
        K_train_nc = K_train;
        K_train_test_nc = K_train_test;
    end

end

