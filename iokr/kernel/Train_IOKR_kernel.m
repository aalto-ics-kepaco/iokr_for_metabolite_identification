function [ train_model ] = Train_IOKR_kernel( KX_list_train, Y_train, ky_param, select_param, iokr_param )
%======================================================
% DESCRIPTION:
% Training of IOKR in the case of a kernel represention in output
%
% INPUTS:
% KX_train:     input kernel matrix on the training set (size n_train*n_train)
% lambda:       regularization parameter (>0)
%
% OUTPUTS:
% C:            (partial) regression model
%
%======================================================
    
    % Selection of the regularization parameter and of the output kernel parameter(s)
    [lambda_opt, KY_par_opt, w_opt] = Select_param_IOKR_kernel2(KX_list_train, Y_train, ky_param, select_param, iokr_param);
    
    % Kernels processing and kernel combination
    [KX_train, process_input] = input_kernel_preprocessing_train(KX_list_train, w_opt, iokr_param);
    
    % Training IOKR with the selected parameter
    n_train = size(KX_train,1);

    C = lambda_opt*eye(n_train) + KX_train;
    
    train_model = struct('C', C, 'KY_par', KY_par_opt, 'w', w_opt, 'process_input', process_input);
    
end

