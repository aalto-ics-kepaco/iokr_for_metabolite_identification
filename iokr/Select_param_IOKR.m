function [ lambda_opt, KY_par_opt, w_opt ] = Select_param_IOKR (KX_list_train, ...
    Y_train, ky_param, opt_param, iokr_param, verbose)
%======================================================
% DESCRIPTION:
% Selection of the regularization parameter in IOKR in the case of a kernel represention in output
%
% INPUTS:
% KX_list_train:    cell array containing the training input kernel matrices
% Y_train:          matrix of size d*n_train containing the training output vectors
% output_param:     1*1 struct array containing the information related to the outputs
% select_param:     1*1 struct array containing information related to the parameter selection 
% iokr_param:       1*1 struct array containing information relative to
%                   centering and multiple kernel learning
%
% OUTPUT:
% lambda_opt:   selected regularization parameter
% KY_par_opt:   struct array containing the selected parameter(s) for the output kernel
% w_opt:        MKL weights when using the output kernel with selected parameter(s)
%
%======================================================    

    if nargin < 6
        verbose = false;
    end % if

    % Possible values for the regularization parameter
    val_lambda = opt_param.val_lambda;
    
    if strcmp(opt_param.cv_type,'cv') && ~isfield(opt_param,'cv_partition')
        % Create cross-validation partition if it is missing
        n_folds = opt_param.num_folds;
        opt_param.cv_partition = cvpartition(size(Y_train,2), 'k', n_folds);
        opt_param.num_folds = n_folds;
    end
    
    switch ky_param.representation
        case 'feature'            
            % Multiple kernel learning
            switch ky_param.type
                case 'linear'
                    KY_train = build_kernel (Y_train, Y_train, ky_param);
                case 'gaussian'
                    assert (strcmp (ky_param.param_selection, 'entropy'), ...
                        'Random fourier features currently only implemented for entropy-gamma-selection.');
                    ky_param.gamma = select_gamma_entropy (Y_train, ky_param);
                    
                    % RANDOM FOURIER FEATURES
                    ky_param.rff = RandomFourierFeatures (size (Y_train, 1), ky_param.rff_dimension);
                    Y_train      = ky_param.rff.getRandomFourierFeatures (Y_train, ky_param.gamma);
                    KY_train     = build_kernel (Y_train, Y_train, struct ('type', 'linear'));
                otherwise
                    error ('Select_param_IOKR:InvalidArgument', ...
                        'Using features as output representation currently only "linear" and "gaussian" kernels are supported. Not %s.', ...
                        ky_param.type);
            end % switch    
            w = mkl_weight(iokr_param.mkl, KX_list_train, normmat (KY_train));
            
            % Input kernels processing and combination
            KX_train = input_kernel_preprocessing_train(KX_list_train, w, iokr_param.center);
            
            % Output feature vectors processing
            Psi_train = output_feature_preprocessing_train(Y_train, iokr_param.center);
            
            % Computation of the MSE for the different regularization parameters
            mse = IOKR_feature_eval_mse(KX_train, Psi_train, opt_param);
            
            % Parameter selection
            [~, ind_lambda_opt] = min(mse);
            lambda_opt = val_lambda(ind_lambda_opt);
            w_opt = w;
            KY_par_opt = ky_param;
            
        case 'kernel'
        
            % If the param_selection field is set to 'entropy', we use the gamma parameter 
            % that maximizes the entropy of the kernel (can only used with Gaussian kernels)
            if isfield(ky_param,'param_selection') && strcmp(ky_param.param_selection,'entropy')
                ky_param.gamma = select_gamma_entropy(Y_train, ky_param);
            end

            % Generates all possible parameter combinations for the output kernel
            ky_param_all_comb = IterGrid(ky_param);

            mse = zeros(length(ky_param_all_comb), length(val_lambda));
            w = cell(length(ky_param_all_comb), 1);
            for ip = 1:length(ky_param_all_comb)
                % Multiple kernel learning
                KY_train = build_kernel(Y_train, Y_train, ky_param_all_comb(ip));
                w{ip} = mkl_weight(iokr_param.mkl, KX_list_train, normmat(KY_train));

                % Input kernels processing and combination
                KX_train = input_kernel_preprocessing_train(KX_list_train, w{ip}, iokr_param.center);

                % Output kernel processing
                KY_train = output_kernel_preprocessing_train(Y_train, ky_param_all_comb(ip), iokr_param.center);

                % Computation of the MSE for the different regularization parameters
                mse(ip,:) = IOKR_kernel_eval_mse(KX_train, KY_train, opt_param, verbose);
            end

            % Parameter selection
            [~, I] = min(mse(:));
            [ind_KY_param_opt, ind_lambda_opt] = ind2sub(size(mse),I);

            lambda_opt = val_lambda(ind_lambda_opt);
            KY_par_opt = ky_param_all_comb(ind_KY_param_opt);
            w_opt = w{ind_KY_param_opt};
    end
end


