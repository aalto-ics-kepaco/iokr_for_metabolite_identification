function [ mse ] = IOKR_feature_eval_mse( KX_train, Psi_train, select_param )
%======================================================
% DESCRIPTION:
% Computation of the mean squared errors (mse) for different values of the 
% regularization parameter
%
% INPUTS:
% KX_train:     training input Gram matrix
% Psi_train:    training output feature vectors
% select_param: 1*1 struct array containing information related to the parameter selection
%
% OUTPUT:
% mse:          vector containing the mse obtained for each regularization parameter
%               in select_param.lambda
%
%====================================================== 

    % Vector containing the possible values for the regularization parameter
    val_lambda = select_param.lambda;

    switch select_param.cv_type
        
        % Parameter selection using inner cross-validation
        case 'cv'
        
            c = select_param.cv_partition; % CV partition
            n_folds = select_param.num_folds; % number of folds

            mse_cv = zeros(length(val_lambda), n_folds);

            for j = 1:n_folds % Cross-validation
                train_set_cv = find(training(c,j));
                test_set_cv = find(test(c,j));
                
                n_train_cv = length(train_set_cv);
                n_test_cv = length(test_set_cv);

                KX_train_cv = KX_train(train_set_cv, train_set_cv);
                KX_train_test_cv = KX_train(train_set_cv, test_set_cv);
                
                Psi_train_cv = Psi_train(:, train_set_cv);
                Psi_test_cv = Psi_train(:, test_set_cv);

                for il = 1:length(val_lambda)

                    % Training                    
                    C_cv = Psi_train_cv / (val_lambda(il)*eye(n_train_cv) + KX_train_cv);
               
                    % Prediction
                    Psi_pred = C_cv * KX_train_test_cv;

                    % Computation of the mean squared error
                    mse_cv(il,j) = 1/n_test_cv*norm(Psi_pred - Psi_test_cv,'fro')^2;
                end
            end
            mse = mean(mse_cv, 2); % mean over the different folds
            
     
        % Parameter selection using leave-one-out cross-validation
        case 'loocv'
        
            n_train = size(Psi_train,2);

            mse = zeros(length(val_lambda),1);
            for il = 1:length(val_lambda)

                B = (val_lambda(il)*eye(n_train) + KX_train) \ KX_train;

                % Computation of the mean squared error
                LOOE = (eye(n_train)-B) / diag(diag(eye(n_train)-B));
                mse(il) = 1/n_train * norm(Psi_train * LOOE,'fro')^2;
            end
    end

end

