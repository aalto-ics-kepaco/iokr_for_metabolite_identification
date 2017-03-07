function [ lambda_opt, KY_par_opt, w_opt ] = Select_param_IOKR_kernel( KX_list_train, Y_train, KY_opt, val_lambda, param )
%======================================================
% DESCRIPTION:
% Selection of the regularization parameter in IOKR in the case of a kernel represention in output
%
% INPUTS:
% KX_list_train:    cell array containing the training input kernel matrices
% Y_train:          matrix of size d*n_train containing the training output vectors
% KY_type:          string indicating the output kernel type ('linear','gaussian' or 'polynomial')
% val_lambda:       vector of strictly positive values among which the regularization 
%                   parameter will be selected
% param:            structure containing the MP-IOKR parameters
%   param.center:   binary value indicating if the input and output
%                   kernel/feature vectors should be centered (1) or not (0)
%   param.mkl:      string indicating the MKL algorithm used for kernel combination 
%                   ('alignf' or 'unimkl')
%   param.cv:       string indicating the type of cross-validation
%                   ('cv' or 'loocv') for parameter selection
%
% OUTPUT:
% lambda_opt:   selected regularization parameter
% KY_par_opt:   selected parameter for the output kernel
% w_opt:        MKL weights learned for the output kernel with parameter KY_par_opt
%
%======================================================    

    % Parameter(s) of the output kernel
    switch KY_opt.type
        case 'gaussian'
            val_KY_param = [0.001 0.005 0.01 0.05 0.1 1];
        case 'polynomial'
            val_KY_param1 = [0 1e-3 1e-2 1e-1 1];
            val_KY_param2 = [1 2 4 6 8];
            [A,B] = meshgrid(val_KY_param1,val_KY_param2);
            val_KY_param = [A(:),B(:)]';
    end

    n_param = size(val_KY_param,2);
    n_train = size(Y_train,2);

    switch param.cv 
        case 'cv'
            n_folds = 10;
            c = cvpartition(n_train, 'k', n_folds); % splits the training examples in n_folds folds

            mse = zeros(n_param, length(val_lamda), n_folds);
            
            for ip = 1:n_param
                
                KY_par = set_kernel_param( KY_opt, val_KY_param(:,ip));
                KY_train = build_kernel(Y_train, Y_train, KY_par);
                
                for j = 1:n_folds
                    train_set_cv = find(training(c,j));
                    test_set_cv = find(test(c,j));
                    
                    % MKL
                    [KX_train_cv, KX_train_test_cv, ~] = mkl(KX_list_train, ...
                        normmat(KY_train(train_set_cv,train_set_cv)), train_set_cv, test_set_cv, param);
                    
                    % Centering and normalization
                    [KY_train_cv, KY_train_test_cv] = input_kernel_center_norm( KY_train, train_set_cv, test_set_cv, param.center);

                    for il = 1:length(val_lambda)
                        lambda = val_lambda(il);

                        % Training and prediction
                        C_cv = Train_IOKR_kernel(KX_train_cv, lambda);
                        B_cv = C_cv \ KX_train_test_cv;


                        % Compute the mean squared error
                        mse(ip,il,j) = 1 + 1/length(test_set_cv)*...
                               trace(B_cv'*KY_train_cv*B_cv - 2*B_cv'*KY_train_test_cv);
                    end
                end
            end
            
            mse_tot = mean(mse,3); % mean over the different folds
            
            [A1,I1] = min(mse_tot,[],1);
            [~,I2] = min(A1,[],2);

            ind_lambda_opt = I2;
            ind_KY_param_opt = I1(:,ind_lambda_opt);

            lambda_opt = val_lambda(ind_lambda_opt);
            
            KY_par_opt = set_kernel_param( KY_opt, val_KY_param(:,ind_KY_param_opt));

            
        case 'loocv'
            
            mse = zeros(n_param, length(val_lambda));
            w = cell(n_param, 1);
            
            for ip = 1:n_param
                
                % MKL
                KY_par = set_kernel_param(KY_opt, val_KY_param(:,ip));
                KY_train = build_kernel(Y_train, Y_train, KY_par);
                
                w{ip} = mkl_weight(param.mkl, KX_list_train, normmat(KY_train));
                
                KX_train = mkl_combine_train(KX_list_train, w{ip}, par.center);
                                
                KY_train_c = center(KY_train, mean(KY_train,1), par.center);
                KY_train_c = normmat(KY_train_c);
            
                for il = 1:length(val_lambda)

                    B = (val_lambda(il)*eye(n_train) + KX_train) \ KX_train;
                    LOOE = (eye(n_train)-B) / diag(diag(eye(n_train)-B));

                    % Compute the mean squared error
                    mse(ip,il) = 1/n_train * trace(LOOE' * KY_train_c * LOOE);
                    
                end
            end
            
        [A1,I1] = min(mse,[],1);
        [~,I2] = min(A1,[],2);

        ind_lambda_opt = I2;
        ind_KY_param_opt = I1(:,ind_lambda_opt);
        
        lambda_opt = val_lambda(ind_lambda_opt);
        KY_param_opt = val_KY_param(:,ind_KY_param_opt);
        KY_par_opt = set_kernel_param(KY_opt, KY_param_opt);
        w_opt = w{ind_KY_param_opt};
    end

end

function [] = selec_cv_without_kernel_param(KX_list_train, KY_train, c, param, val_lambda)


    for j = 1:n_folds
        train_set_cv = find(training(c,j));
        test_set_cv = find(test(c,j));

        % MKL
        [KX_train_cv, KX_train_test_cv, ~] = mkl(KX_list_train, ...
            normmat(KY_train(train_set_cv,train_set_cv)), train_set_cv, test_set_cv, param);

        % Centering and normalization
        [KY_train_cv, KY_train_test_cv] = input_kernel_center_norm( KY_train, train_set_cv, test_set_cv, param.center);

        for il = 1:length(val_lambda)
            lambda = val_lambda(il);

            % Training and prediction
            C_cv = Train_IOKR_kernel(KX_train_cv, lambda);
            B_cv = C_cv \ KX_train_test_cv;


            % Compute the mean squared error
            mse(ip,il,j) = 1 + 1/length(test_set_cv)*...
                   trace(B_cv'*KY_train_cv*B_cv - 2*B_cv'*KY_train_test_cv);
        end
    end
end
