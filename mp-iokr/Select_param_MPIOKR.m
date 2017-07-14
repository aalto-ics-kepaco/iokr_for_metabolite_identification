function [ lambda_opt, gamma_opt, KY_par_opt, w_opt ] = Select_param_MPIOKR ( ...
    KX_list_train, Y_train, Y_C_train, ky_param, opt_param, mp_iokr_param )
%======================================================
% DESCRIPTION:
% Selection of the regularization parameter in MP-IOKR in the case of a kernel represention in output
%
% INPUTS:
% KX_list_train:    cell array containing the training input kernel matrices
% Y_train:          matrix of size d*n_train containing the training output vectors
% Y_C_train:        cell array containing the training candidate output vectors
% ky_param:         1*1 struct array containing the information related to the output kernel
% select_param:     1*1 struct array containing information related to the parameter selection 
% mp_iokr_param:    1*1 struct array containing information relative to
%                   centering and multiple kernel learning
%
% OUTPUT:
% lambda_opt, gamma_opt: selected regularization parameters
% KY_par_opt:       struct array containing the selected parameter(s) for the output kernel
% w_opt:            MKL weights when using the output kernel with selected parameter(s)
%
%======================================================    

    % Possible values for the regularization parameter
    val_lambda = opt_param.val_lambda;
    val_gamma  = opt_param.val_gamma;
    
    n_kx = length(KX_list_train); % number of input kernels
    
    % Create cross-validation partition
    opt_param.cv = cvpartition (size(Y_train,2), 'k', opt_param.nInnerFolds);

    % If the param_selection field is set to 'entropy', we use the gamma parameter 
    % that maximizes the entropy of the kernel (can only used with Gaussian kernels)
    if strcmp (ky_param.type, 'gaussian')
        if isfield (ky_param, 'param_selection') && strcmp (ky_param.param_selection, 'entropy')
            ky_param.gamma = select_gamma_entropy (Y_train, ky_param);
        end
    end
         
    % Generates all possible parameter combinations for the output kernel
    ky_param_all_comb = IterGrid(ky_param);

    %% Parameter selection    
    mse = zeros (length (val_lambda), length (ky_param_all_comb));
    w = cell (length (ky_param_all_comb), 1);
    for ip = 1:length(ky_param_all_comb)
        % Multiple kernel learning
        KY_train = build_kernel (Y_train, Y_train, ky_param_all_comb(ip));
        w{ip} = mkl_weight (mp_iokr_param.mkl, KX_list_train, normmat (KY_train));

        % Input kernels processing and combination
        KX_train = mpiokr_input_kernel_preprocessing_train (KX_list_train, w{ip}, mp_iokr_param);

        % Output processing
        if strcmp (ky_param.type, 'linear')
            Psi_train = output_feature_preprocessing_train (Y_train, mp_iokr_param.center);
        else
            KY_train = output_kernel_preprocessing_train (Y_train, ky_param_all_comb(ip), mp_iokr_param.center);
        end

        % Selection of the regularization parameter of reverse IOKR
        gamma_opt = zeros(n_kx,1);
        for i = 1:n_kx
            if strcmp(ky_param.type,'linear')
                mse_gamma = IOKR_feature_eval_mse(KX_train{i}, Psi_train, struct('lambda', val_gamma, 'cv_type', 'loocv'));
            else
                mse_gamma = IOKR_kernel_eval_mse(KX_train{i}, KY_train, struct('lambda', val_gamma, 'cv_type', 'loocv'));
            end
            [~,ind_gamma_opt] = min(mse_gamma);
            gamma_opt(i) = val_gamma(ind_gamma_opt);      
        end

        % Modify the KX_list in order to combine the kernel belonging to unique gamma values
        if (strcmp (mp_iokr_param.rev_iokr, 'separate'))
            KX_train = mpiokr_input_kernel_preprocessing_train (KX_list_train, w, mp_iokr_param, gamma_opt);
            gamma_opt_red = unique (gamma_opt);
            n_kx = length (gamma_opt_red);
        else
            gamma_opt_red = gamma_opt;
        end

        Mc = cell(n_kx,1);
        for k = 1:n_kx
            if strcmp(ky_param.type,'linear')
                Mc{k} = (gamma_opt_red(k) * eye(n_train) + (Psi_train'*Psi_train)) \ (Psi_train');
            else
                % TODO: Can we use the Cholesky decomposition here?
                Mc{k} = inv (gamma_opt_red(k) * eye(n_train) + KY_train);
            end
        end
        M = cell2mat (Mc);

        % Computation of the MSE for the different regularization parameters
        if strcmp (ky_param.type, 'linear')
            mse(:,ip) = MPIOKR_feature_eval_mse (KX_train, Psi_train, Y_C_train, M, mp_iokr_param, opt_param);
        else
            mse(:,ip) = MPIOKR_kernel_eval_mse (KX_train, Y_train, Y_C_train, ky_param_all_comb(ip), M, mp_iokr_param, opt_param);
        end
    end

    % Parameter selection
    [~, I] = min(mse(:));
    [ind_KY_param_opt, ind_lambda_opt] = ind2sub(size(mse),I);

    lambda_opt = val_lambda(ind_lambda_opt);
    KY_par_opt = ky_param_all_comb(ind_KY_param_opt);
    w_opt = w{ind_KY_param_opt};

end

