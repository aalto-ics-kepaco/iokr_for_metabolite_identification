function [ train_model ] = Train_IOKR( KX_list_train, Y_train, ...
    ky_param, opt_param, iokr_param, verbose )
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
    if nargin < 6
        verbose = false;
    end % if
    
    % Selection of the regularization parameter and of the output 
    % kernel parameter(s)
    if verbose; tic; end % if
    [lambda_opt, KY_par_opt, w_opt] = Select_param_IOKR ( ...
        KX_list_train, Y_train, ky_param, opt_param, iokr_param, verbose);
    if verbose 
        fprintf('IOKR parameter selection took %.3fs.\n', toc);
        fprintf('Optimal regularization parameter: %f.\n', lambda_opt);
    end % if

    % Kernels processing and kernel combination
    if verbose; tic; end % if
    [KX_train, process_input] = input_kernel_preprocessing_train (...
        KX_list_train, w_opt, iokr_param.center);
    if verbose 
        fprintf('Full data input kernel processing took %.3fs.\n', toc);
    end % if
    
    % Training IOKR with the selected parameters
    if verbose; tic; end % if
    switch iokr_param.model_representation
        case 'only_C'
            C = lambda_opt*eye(size(KX_train)) + KX_train;
        case 'Chol_decomp_of_C'
            C = chol (lambda_opt*eye(size(KX_train)) + KX_train, 'lower');
    end % switch     
    if verbose
        fprintf('Calculation of the model using all training data took %.3fs.\n', toc);
    end % if
    
    train_model = struct (                               ...
        'C',                    C,                       ...
        'ker_center',           iokr_param.center,       ...
        'lambda_opt',           lambda_opt,              ...
        'process_input',        process_input,           ...
        'KY_par',               KY_par_opt,              ...
        'representation',       ky_param.representation, ...
        'model_representation', iokr_param.model_representation);
end

