function [ score ] = IOKR_kernel(KX_list, train_set, test_set, Y_train, Y_C_test, val_lambda, param )
%======================================================
% DESCRIPTION:
% IOKR in the case of a kernel represention in output
%
% INPUTS:
% KX_list:      list of input kernel matrices (each of size n*n)
% train_set:    vector containing the training indices (length n_train)
% test_test:    vector containing the test indices (length n_test)
% Y_train:      matrix of size d*n_train containing the training output vectors
% Y_C_test:     test candidate sets
% val_lambda:   vector of strictly positive values among which the regularization 
%               parameter of IOKR will be selected
% param:        structure containing the MP-IOKR parameters
%   param.center:       binary value indicating if the input and output
%                       kernel/feature vectors should be centered (1) or not (0)
%   param.mkl:          string indicating the MKL algorithm used for kernel combination 
%                       ('alignf' or 'unimkl')
%   param.cv:           string indicating the type of cross-validation
%                       ('cv' or 'loocv') for parameter selection
%
% OUTPUTS:
% score:        cell of length n_test, each entry being a vector containing 
%               the scores associated with the candidates in the corresponding
%               candidate set
%
%======================================================

    KX_list_train = cellfun(@(x) x(train_set,train_set), KX_list, 'UniformOutput', false);
    
    % Selection of the regularization parameter and of the output kernel parameter(s)
    [lambda_opt, KY_par_opt, w_opt] = Select_param_IOKR_kernel(KX_list_train, Y_train, KY_type, val_lambda, param);
    
    % Input kernel combination
    [KX_train, KX_train_test] = mkl_combine_train_test(KX_list, train_set, test_set, w_opt, param.center);

    % Training IOKR with the selected parameter
    C = Train_IOKR_kernel(KX_train, lambda_opt);
    B = C \ KX_train_test;
    
    % Preimage
    score = Preimage_IOKR_kernel(Y_train, KY_par_opt, B, Y_C_test, param);

end

