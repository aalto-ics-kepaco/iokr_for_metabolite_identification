function [ train_model ] = Train_IOKR( KX_list_train, Y_train, output_param, select_param, iokr_param )
%======================================================
% DESCRIPTION:
% Training step of IOKR
%
% INPUTS:
% KX_list_train:    cell array containing the training input kernel matrices
% Y_train:          matrix of size d*n_train containing the training output vectors
% output_param:     1*1 struct array containing the information related to the outputs
% select_param:     1*1 struct array containing information related to the parameter selection
% iokr_param:       1*1 struct array containing information relative to
%                   centering and multiple kernel learning
%
% OUTPUTS:
% train_model:      1*1 struct array containing the regression model
%                   and information on the training data 
%
%======================================================
    
    % Selection of the regularization parameter and of the output kernel parameter(s)
    [lambda_opt, KY_par_opt, w_opt] = Select_param_IOKR(KX_list_train, Y_train, output_param, select_param, iokr_param);
    
    % Kernels processing and kernel combination
    [KX_train, process_input] = input_kernel_preprocessing_train(KX_list_train, w_opt, iokr_param.center);
    
    % Training IOKR with the selected parameter
    C = lambda_opt*eye(size(KX_train)) + KX_train;
    
    train_model = struct('C', C, 'process_input', process_input, 'KY_par', KY_par_opt, ...
                    'representation', output_param.representation);
    
end

