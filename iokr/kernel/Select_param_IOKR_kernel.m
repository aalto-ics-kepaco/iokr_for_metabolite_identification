function [ lambda_opt, KY_par_opt, w_opt ] = Select_param_IOKR_kernel( KX_list_train, Y_train, ky_param, select_param, iokr_param )
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

    % If the param_selection field is set to 'entropy', we use the gamma parameter 
    % that maximizes the entropy of the kernel (can only used with Gaussian kernels)
    if isfield(ky_param,'param_selection') && strcmp(ky_param.gamma_selection,'entropy')
        ky_param.gamma = select_gamma_entropy(Y_train, ky_param);
    end
    
    if strcmp(select_param.cv_type,'cv') && ~isfield(select_param,'cv_partition')
        % Create cross-validation partition
        n_folds = select_param.num_folds;
        select_param.cv_partition = cvpartition(size(Y_train,2), 'k', n_folds);
    end

    % Possible values for the regularization parameter
    val_lambda = select_param.lambda;
    
    % Generates all possible parameter combinations for the output kernel
    ky_param_all_comb = IterGrid( ky_param );
    
    mse = zeros(length(ky_param_all_comb), length(val_lambda));
    w = cell(length(ky_param_all_comb), 1);
    
    for ip = 1:length(ky_param_all_comb)

        KY_train = build_kernel(Y_train, Y_train, ky_param_all_comb(ip));
        
        % Multiple kernel learning and kernel centering/normalization
        w{ip} = mkl_weight(iokr_param.mkl, KX_list_train, normmat(KY_train));
        KX_train = mkl_combine_train(KX_list_train, w{ip}, iokr_param.center);
        
        % Centering and normalization
        KY_train_c = center(KY_train, mean(KY_train,1), iokr_param.center);
        KY_train_c = normmat(KY_train_c);

        % Computation of the mse for different regularization parameters lambda
        mse(ip,:) = IOKR_eval_mse(KX_train, KY_train_c, val_lambda, select_param);
    end
       
    [A1,I1] = min(mse,[],1);
    [~,I2] = min(A1,[],2);
    
    ind_lambda_opt = I2;
    ind_KY_param_opt = I1(:,ind_lambda_opt);
    
    lambda_opt = val_lambda(ind_lambda_opt);
    KY_par_opt = ky_param_all_comb(ind_KY_param_opt);
    w_opt = w{ind_KY_param_opt};

end


function [ mse ] = IOKR_eval_mse( KX_train, KY_train, val_lambda, select_param )

    switch select_param.cv_type
        
        % Parameter selection with inner cross-validation
        case 'cv'
        
            c = select_param.cv_partition;
            n_folds = c.NumTestSets;

            mse_cv = zeros(1, length(val_lambda), n_folds);

            for j = 1:n_folds
                train_set_cv = find(training(c,j));
                test_set_cv = find(test(c,j));

                KX_train_cv = KX_train(train_set_cv, train_set_cv);
                KX_train_test_cv = KX_train(train_set_cv, test_set_cv);

                KY_train_cv = KY_train(train_set_cv, train_set_cv);
                KY_train_test_cv = KY_train(train_set_cv, test_set_cv);

                for il = 1:length(val_lambda)

                    % Training and prediction
                    C_cv = Train_IOKR_kernel(KX_train_cv, val_lambda(il));

                    B_cv = Prediction_IOKR_kernel(KX_train_test_cv, C_cv);

                    % Compute the mean squared error
                    mse_cv(1,il,j) = 1 + 1/length(test_set_cv)*...
                           trace(B_cv'*KY_train_cv*B_cv - 2*B_cv'*KY_train_test_cv);
                end
            end
            mse = mean(mse_cv, 3); % mean over the different folds
     
        % Parameter selection with leave-one-out cross-validation
        case 'loocv'
        
            n_train = size(KY_train,1);

            mse = zeros(1, length(val_lambda));
            for il = 1:length(val_lambda)

                % Training
                C = Train_IOKR_kernel(KX_train, val_lambda(il));

                % Prediction on the training set
                B = Prediction_IOKR_kernel(KX_train, C);

                % Compute the mean squared error
                LOOE = (eye(n_train)-B) / diag(diag(eye(n_train)-B));
                mse(1,il) = 1/n_train * trace(LOOE' * KY_train * LOOE);
            end
    end

end
