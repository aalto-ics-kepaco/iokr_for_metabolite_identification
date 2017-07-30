function [ vec_gamma_opt ] = Select_param_reverse_IOKR( KX_train_list, Psi_train, opt_param )
%======================================================
% DESCRIPTION:
% Selection of the regularization parameters for the reverse IOKR approach 
% in the case of multiple input kernels using leave-one-out cross-validation
% (one parameter per input kernel is selected)
%
% INPUTS:
% KX_train_list:    list of input Gram matrices (size n_train*n_train)
% Psi_train:        matrix of size d*n_train containing the training output feature vectors
% val_gamma:        vector of stricly positive values among which the regularization 
%                   parameters of reverse IOKR will be selected
%
% OUTPUTS:
% vec_gamma_opt:    vector containing the different selected
%                   regularization parameters (one for each input kernel)
%
%======================================================

    n_kx = length(KX_train_list); % number of input kernels
        
    vec_gamma_opt = zeros(n_kx,1);
    for i = 1:n_kx
        
        mse = IOKR_kernel_eval_mse(Psi_train'*Psi_train, KX_train_list{i}, ...
            struct('val_lambda', opt_param.val_gamma, 'cv_type', 'loocv'));
        
        [~,ind_gamma_opt] = min(mse);
        vec_gamma_opt(i) = opt_param.val_gamma(ind_gamma_opt);      
    end
end

