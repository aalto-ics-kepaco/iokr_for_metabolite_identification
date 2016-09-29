function [ score, debug_info ] = IOKR_feat( KX_list, Y_train, Y_C_test, ...
    opt_param, iokr_param, data_param, debug_param)
%======================================================
% DESCRIPTION:
% IOKR in the case of a feature represention in output
%
% INPUTS:
% KX_list:      list of input kernel matrices (each of size n*n)
% train_set:    vector containing the training indices (length n_train)
% test_test:    vector containing the test indices (length n_test)
% Y_train:      matrix of size d*n_train containing the training output feature vectors
% Y_C_test:     test candidate sets
% opt_param:    structure containing the optimization parameters
%
%   opt_param.val_gamma:  vector of strictly positive values among which 
%                         the regularization parameter of reverse IOKR will 
%                         be selected
% iokr_param:   structure containing the MP-IOKR parameters
%   iokr_param.center:  binary value indicating if the input and output
%                       kernel/feature vectors should be centered (1) or not (0)
%   iokr_param.mkl:     string indicating the MKL algorithm used for kernel combination 
%                       ('alignf' or 'unimkl')
%   iokr_param.cv_type: string indicating the type of cross-validation
%                       ('cv' or 'loocv') for parameter selection
% data_param:   structure containing the parameters for the data handling
%
%   data_param.train_set        binary vector containing indicating the
%                               examples used for training (length n, 
%                               sum (data_param.train_set) == n_train)
%   data_param.test_set         binary vector containing indicating the
%                               examples used for testing (length n, 
%                               sum (data_param.test_set) == n_test)
%
% OUTPUTS:
% score:        cell of length n_test, each entry being a vector containing 
%               the scores associated with the candidates in the corresponding
%               candidate set
%
%======================================================
    train_set = find (data_param.train_set);
    test_set  = find (data_param.test_set);
    
    % Learning kernel weights with Multiple Kernel Learning  
    KX_list_train = cellfun(@(x) x(train_set,train_set), KX_list, 'UniformOutput', false);
    w = mkl_weight(iokr_param.mkl, KX_list_train, normmat(Y_train'*Y_train));
    clear KX_list_train
    
    % Combining kernels
    [KX_train, KX_train_test] = mkl_combine_train_test(KX_list, train_set, test_set, w, iokr_param.center);
    
    % Training output feature vectors
    mean_Y_train = mean(Y_train,2);
    Psi_train = norma(Y_train, mean_Y_train, iokr_param.center);
    
    %% Parameter selection
    fprintf ('Parameter selection\n');
    debug_info = struct ();
    
    [lambda_opt, debug_info.err] = Select_param_IOKR_feat(KX_train, Psi_train, ...
        opt_param.val_lambda, iokr_param.cv_type);
    
    debug_info.lambda_opt = lambda_opt;
    
    fprintf ('Selected parameter: lambda = %e\n', debug_info.lambda_opt);
    
    if (debug_param.verbose)
        fprintf ('MSE:\n')
        disp (debug_info.err);
    end % if
    
    % Training
    C = Train_IOKR_feat(KX_train, Psi_train, lambda_opt);
    
    % Prediction on the test set
    Psi_pred = Prediction_IOKR_feat(KX_train_test, C);
    
    % Preimage
    score = Preimage_IOKR_feat(Psi_pred, Y_C_test, mean_Y_train, iokr_param.center);
end

