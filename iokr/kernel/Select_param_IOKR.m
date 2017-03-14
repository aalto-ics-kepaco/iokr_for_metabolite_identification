function [ lambda_opt, KY_par_opt, w_opt ] = Select_param_IOKR( KX_list_train, Y_train, output_param, select_param, iokr_param )
%======================================================
% DESCRIPTION:
% Selection of the regularization parameter in IOKR in the case of a kernel represention in output
%
% INPUTS:
% KX_list_train:    cell array containing the training input kernel matrices
% Y_train:          matrix of size d*n_train containing the training output vectors
% ky_param:         1*1 struct array containing the information related to the
%                   output kernel
% select_param:     1*1 struct array containing information related to the parameter selection 
% iokr_param:            structure containing the IOKR parameters
%   param.center:   binary value indicating if the input and output
%                   kernel/feature vectors should be centered (1) or not (0)
%   param.mkl:      string indicating the MKL algorithm used for kernel combination 
%                   ('alignf' or 'unimkl')
%   param.cv_type:       string indicating the type of cross-validation
%                   ('cv' or 'loocv') for parameter selection
%
% OUTPUT:
% lambda_opt:   selected regularization parameter
% KY_par_opt:   struct array containing the selected parameter(s) for the output kernel
% w_opt:        MKL weights when using the output kernel with selected
%               parameter(s)
%
%======================================================    

    % Possible values for the regularization parameter
    val_lambda = select_param.lambda;
    
    if strcmp(select_param.cv_type,'cv') && ~isfield(select_param,'cv_partition')
        % Create cross-validation partition if it is missing
        n_folds = select_param.num_folds;
        select_param.cv_partition = cvpartition(size(Y_train,2), 'k', n_folds);
        select_param.num_folds = n_folds;
    end
    
    switch output_param.representation
        case 'feature'
                        
            % Multiple kernel learning
            KY_par = struct('type','linear');
            KY_train = build_kernel(Y_train, Y_train, KY_par);
            w = mkl_weight(iokr_param.mkl, KX_list_train, normmat(KY_train));
            
            % Input kernels processing and combination
            KX_train = input_kernel_preprocessing_train(KX_list_train, w, iokr_param);
            
            % Output feature vectors processing
            Psi_train = output_feature_preprocessing_train(Y_train, iokr_param.center);
            
            % Computation of the MSE for the different regularization parameters
            mse = IOKR_feature_eval_mse(KX_train, Psi_train, select_param);
            
            % Parameter selection
            [~, ind_lambda_opt] = min(mse);
            lambda_opt = val_lambda(ind_lambda_opt);
            w_opt = w;
            KY_par_opt = KY_par;
            
        case 'kernel'
            
            ky_param = output_param.kernel_param;

            % If the param_selection field is set to 'entropy', we use the gamma parameter 
            % that maximizes the entropy of the kernel (can only used with Gaussian kernels)
            if isfield(ky_param,'param_selection') && strcmp(ky_param.gamma_selection,'entropy')
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
                mse(ip,:) = IOKR_kernel_eval_mse(KX_train, KY_train, select_param);
            end

            % Parameter selection
            [~, I] = min(mse(:));
            [ind_KY_param_opt, ind_lambda_opt] = ind2sub(size(mse),I);

            lambda_opt = val_lambda(ind_lambda_opt);
            KY_par_opt = ky_param_all_comb(ind_KY_param_opt);
            w_opt = w{ind_KY_param_opt};
    end

end


