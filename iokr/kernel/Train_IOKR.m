function [ train_model ] = Train_IOKR( KX_list_train, Y_train, output_param, select_param, iokr_param )
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
    [lambda_opt, KY_par_opt, w_opt] = Select_param_IOKR_kernel(KX_list_train, Y_train, output_param, select_param, iokr_param);
    
    % Kernels processing and kernel combination
    [KX_train, process_input] = input_kernel_preprocessing_train(KX_list_train, w_opt, iokr_param.center);
    
    [~, process_output] = output_kernel_preprocessing_train(Y_train, KY_par_opt, iokr_param.center);
    
    % Training IOKR with the selected parameter
    C = lambda_opt*eye(size(KX_train)) + KX_train;
    
    train_model = struct('C', C, 'process_input', process_input, 'process_output', process_output, 'KY_par', KY_par_opt);
    
end

