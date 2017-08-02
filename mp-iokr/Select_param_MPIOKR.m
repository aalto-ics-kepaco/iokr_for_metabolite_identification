function [ lambda_opt, gamma_opt, KY_par_opt, w_opt ] = Select_param_MPIOKR ( ...
    KX_list_train, Y_train, Y_C_train, ky_param, opt_param, mp_iokr_param, debug_param )
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

    n_train = size (Y_train, 2);

    % Possible values for the regularization parameter
    val_lambda = opt_param.val_lambda;
    val_gamma  = opt_param.val_gamma;
    
    n_kx = length(KX_list_train); % number of input kernels
    
    % Create cross-validation partition
    opt_param.cv_partition = cvpartition (n_train, 'k', opt_param.nInnerFolds);

    % If the param_selection field is set to 'entropy', we use the gamma parameter 
    % that maximizes the entropy of the kernel (can only used with Gaussian kernels)
    if (strcmp (ky_param.type, 'gaussian')) && (strcmp (ky_param.param_selection, 'entropy'))
        ky_param.gamma = select_gamma_entropy (Y_train, ky_param);
    end % if 
    if (strcmp (ky_param.type, 'gaussian')) && (strcmp (ky_param.representation, 'feature'))
        % RANDOM FOURIER FEATURES
        assert (isfield (ky_param, 'rff'), 'Random fourier feature matrix must be provided.');
    end % if
         
    % Generates all possible parameter combinations for the output kernel
    ky_param_all_comb = IterGrid(ky_param);

    %% Parameter selection    
    mse = zeros (length (val_lambda), length (ky_param_all_comb));
    w = cell (length (ky_param_all_comb), 1);
    for ip = 1:length(ky_param_all_comb)
        % Multiple kernel learning
        if (strcmp (ky_param.representation, 'feature')) && (strcmp (ky_param.type, 'gaussian'))
            % RANDOM FOURIER FEATURES
            Psi_train = ky_param.rff.getRandomFourierFeatures (Y_train, ky_param_all_comb(ip).gamma);
            
            % We can approximate the gaussian-features using random fourier
            % features. Calculating the _linear kernel_ of a set of random 
            % fourier features approximates the a gaussian kernel.
            KY_train = build_kernel (Psi_train, Psi_train, struct ('type', 'linear'));
        else
            KY_train = build_kernel (Y_train, Y_train, ky_param_all_comb(ip));
        end % if
        
        w{ip} = mkl_weight (mp_iokr_param.mkl, KX_list_train, normmat (KY_train));

        % Input kernels processing and combination
        KX_train = mpiokr_input_kernel_preprocessing_train (KX_list_train, w{ip}, mp_iokr_param);

        % Output processing
        if strcmp (ky_param.representation, 'feature')
            switch ky_param.type
                case 'linear'
                    [Psi_train, process_output] = output_feature_preprocessing_train ( ...
                        Y_train, mp_iokr_param.center);
                case 'gaussian'
                    % RANDOM FOURIER FEATURES
                    [Psi_train, process_output] = output_feature_preprocessing_train ( ...
                        Psi_train, mp_iokr_param.center);
                otherwise
                    error ('Select_param_MPIOKR:InvalidArgument', ...
                        'Using features as output representation currently only "linear" and "gaussian" kernels are supported. Not %s.', ...
                        ky_param.type)
            end % switch
        else
            KY_train = output_kernel_preprocessing_train (Y_train, ky_param_all_comb(ip), mp_iokr_param.center);
        end

        % Selection of the regularization parameter of reverse IOKR
        gamma_opt = zeros(n_kx,1);
        for i = 1:n_kx
            if strcmp (ky_param.representation, 'feature')
                mse_gamma = IOKR_kernel_eval_mse (build_kernel (Psi_train, Psi_train, struct ('type', 'linear')), ...
                    KX_train{i}, struct('val_lambda', val_gamma, 'cv_type', 'loocv'));
            else
                mse_gamma = IOKR_kernel_eval_mse(KY_train, KX_train{i}, ...
                    struct('val_lambda', val_gamma, 'cv_type', 'loocv'));
            end
            [~,ind_gamma_opt] = min(mse_gamma);
            gamma_opt(i) = val_gamma(ind_gamma_opt);      
        end

        % Modify the KX_list in order to combine the kernel belonging to unique gamma values
        if (strcmp (mp_iokr_param.rev_iokr, 'separate'))
            KX_train = mpiokr_input_kernel_preprocessing_train (KX_list_train, w{ip}, mp_iokr_param, gamma_opt);
            gamma_opt_red = unique (gamma_opt);
            n_kx = length (gamma_opt_red);
        else
            gamma_opt_red = gamma_opt;
        end

        % Computation of the MSE for the different regularization parameters
        if strcmp (ky_param.representation, 'feature')
            mse(:,ip) = MPIOKR_feature_eval_mse (KX_train, Psi_train, process_output.mean, ...
                Y_C_train, gamma_opt_red, mp_iokr_param, opt_param, ky_param_all_comb(ip), debug_param);
        else
            mse(:,ip) = MPIOKR_kernel_eval_mse (KX_train, Y_train, Y_C_train, ...
                ky_param_all_comb(ip), gamma_opt_red, mp_iokr_param, opt_param, debug_param);
        end
    end

    % Parameter selection
    [~, I] = min(mse(:));
    [ind_lambda_opt, ind_KY_param_opt] = ind2sub(size(mse),I);

    lambda_opt = val_lambda(ind_lambda_opt);
    
    if (debug_param.verbose)
        fprintf ('Optimal lambda: %f\n', lambda_opt);
    end % if
    
    KY_par_opt = ky_param_all_comb(ind_KY_param_opt);
    w_opt = w{ind_KY_param_opt};

end

