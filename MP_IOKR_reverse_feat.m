function [ score ] = MP_IOKR_reverse_feat(KX_list, train_set, test_set, Y_train, Y_C, val_gamma, val_lambda, param)
%======================================================
% DESCRIPTION:
% MP-IOKR with reverse IOKR in the case of a feature represention in output
%
% INPUTS:
% KX_list:      list of input kernel matrices (each of size n*n)
% train_set:    vector containing the training indices (length n_train)
% test_test:    vector containing the test indices (length n_test)
% Y_train:      matrix of size d*n_train containing the training output feature vectors
% Y_C:          Candidate sets
% val_gamma:    vector of strictly positive values among which the regularization 
%               parameter of reverse IOKR will be selected
% val_lambda:   vector of strictly positive values among which the regularization 
%               parameter of MP-IOKR will be selected
% param:        structure containing the MP-IOKR parameters
%   param.center:       binary value indicating if the input and output
%                       kernel/feature vectors should be centered (1) or not (0)
%   param.mkl:          string indicating the MKL algorithm used for kernel combination 
%                       ('alignf' or 'unimkl')
%   param.rev_iokr:     string indicating if the input feature maps of the different
%                       input kernels should be learnt jointly ('joint') or separately ('separate')
%
% OUTPUTS:
% score:        cell of length n_test, each entry being a vector containing 
%               the scores associated with the candidates in the corresponding
%               candidate set
%
%======================================================

    ker_center = param.center; % centering option

    % Learning kernel weights with Multiple Kernel Learning  
    KX_list_train = cellfun(@(x) x(train_set,train_set), KX_list, 'UniformOutput', false);
    w = mkl_weight(param.mkl, KX_list_train, normmat(Y_train'*Y_train));
    clear KX_list_train
    
    % Centering and normalization of input kernel
    n_kx = length(KX_list);
    switch param.rev_iokr
        case 'joint' 
            % Computation of the combined input kernel
            [KX_train_combined, KX_train_test_combined] = mkl_combine_train_test(KX_list, train_set, test_set, w, ker_center);
            KX_train_list = {KX_train_combined};
            KX_train_test_list = {KX_train_test_combined};
            clear KX_train_combined KX_train_test_combined
        case 'separate'
            KX_train_list = cell(n_kx,1);
            KX_train_test_list = cell(n_kx,1);
            for k = 1:n_kx
                [KX_train_list{k}, KX_train_test_list{k}] = input_kernel_center_norm(KX_list{k}, train_set, test_set, ker_center);
                KX_train_list{k} = w(k) * KX_train_list{k};
                KX_train_test_list{k} = w(k) * KX_train_test_list{k};
            end
    end

    % Training output feature vectors
    mean_Y_train = mean(Y_train,2);
    Psi_train = norma(Y_train, mean_Y_train, ker_center);

    Y_C_train = Y_C(train_set);
    Y_C_test = Y_C(test_set);

    % Parameter selection
    [gamma_opt, lambda_opt] = Select_param_MP_IOKR_reverse_feat(KX_train_list, Y_train,...
        Y_C_train, val_gamma, val_lambda, ker_center);
    
    % Training the reverse IOKR model
    M = Train_reverse_IOKR(Psi_train, gamma_opt);
    
    % Training the MP-IOKR model
    [Mean_Psi_C_train, Cov_Psi_C_train] = Compute_cov_mean_feat(Y_C_train, mean_Y_train, ker_center);
    C = Train_MP_IOKR_reverse_feat(KX_train_list, Psi_train, M, Mean_Psi_C_train, Cov_Psi_C_train, lambda_opt);
    
    % Prediction on the test set
    Psi_pred = Prediction_MP_IOKR_reverse_feat(KX_train_test_list, C);
    
    % Preimage
    score = Preimage_MP_IOKR_feat(Psi_pred, Y_C_test, mean_Y_train, ker_center);

end

