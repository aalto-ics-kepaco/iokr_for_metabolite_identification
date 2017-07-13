function [ train_model ] = Train_IOKR( KX_list_train, Y_train, ...
    output_param, select_param, iokr_param )
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
    
%     t = cputime; 
    
    % Selection of the regularization parameter and of the output 
    % kernel parameter(s)
    [lambda_opt, KY_par_opt, w_opt] = Select_param_IOKR ( ...
        KX_list_train, Y_train, output_param, select_param, iokr_param);
%     fprintf ('Optimal lambda = %f\n', lambda_opt);
        
%     fprintf ('Parameter selection (CPU-time): %f\n', cputime - t);

    % Kernels processing and kernel combination
    [KX_train, process_input] = input_kernel_preprocessing_train (...
        KX_list_train, w_opt, iokr_param.center);
    
%     t = cputime;
    
    % Training IOKR with the selected parameter
    switch iokr_param.model_representation
        case 'only_C'
            C = lambda_opt*eye(size(KX_train)) + KX_train;
        case 'Chol_decomp_of_C'
            C = chol (lambda_opt*eye(size(KX_train)) + KX_train, 'lower');
    end % switch     
    
%     fprintf ('Training (CPU-time): %f\n', cputime - t);
    
    train_model = struct('C', C, 'process_input', process_input, ...
        'KY_par', KY_par_opt, ...
        'representation', output_param.representation, ...
        'model_representation', iokr_param.model_representation);
end

