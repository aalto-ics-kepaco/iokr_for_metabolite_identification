function [ lambda_opt, mse ] = Select_param_IOKR_feat( KX_train, Psi_train, val_lambda, cv_type )
%======================================================
% DESCRIPTION:
% Selection of the regularization parameter in IOKR
%
% INPUTS:
% KX_train:     input Gram matrix on the training set (size n_train*n_train)
% Psi_train:    matrix containing the output feature vectors (size d*n_train)
% val_lambda:   vector of striclty positive values among which the regularization 
%               parameter will be selected
% cv_type:      string indicating the type of cross-validation ('cv' or 'loocv')
%
% OUTPUT:
% lambda_opt:   selected regularization parameter
%
%======================================================

    if nargin <= 3
        cv_type = 'loocv';
    end
    
    if nargin <= 2
        val_lambda = 10.^(-4:4);
    end
    
    if nargin <= 1
        disp('An argument is missing.')
    end

    n_train = size(Psi_train,2);

    switch cv_type 
        case 'cv'
            n_folds = 10;
            c = cvpartition(n_train, 'k', n_folds); % splits the training examples in n_folds folds

            mse = zeros(length(val_lambda), n_folds);
            for j = 1:n_folds
                
                % Defining training and test sets
                train_set_cv = find(training(c,j));
                test_set_cv = find(test(c,j));
                
                KX_train_cv = KX_train(train_set_cv, train_set_cv);
                KX_train_test_cv = KX_train(train_set_cv, test_set_cv);
                
                Psi_train_cv = Psi_train(:,train_set_cv);
                Psi_test_cv = Psi_train(:,test_set_cv);
                

                for il = 1:length(val_lambda)
                    lambda = val_lambda(il);
                    
                    % Training and prediction
                    C_cv = Train_IOKR_feat(KX_train_cv, Psi_train_cv, lambda);
                    Psi_pred = Pred_IOKR_feat(C_cv, KX_train_test_cv);
                    
                    % Compute the mean squared error
                    mse(il,j) = 1/length(test_set_cv) * norm(Psi_pred - Psi_test_cv,'fro')^2;
                end
            end
            
            mse_tot = mean(mse,2); % mean over the different folds
            [~, ind_lambda_opt] = min(mse_tot);
            lambda_opt = val_lambda(ind_lambda_opt);
            
        case 'loocv'
            mse = zeros(length(val_lambda),1);
            for il = 1:length(val_lambda)

                B = (val_lambda(il)*eye(n_train) + KX_train) \ KX_train;
                LOOE = (eye(n_train)-B) / diag(diag(eye(n_train)-B));

                % Compute the mean squared error
                mse(il) = 1/n_train * norm(Psi_train * LOOE,'fro')^2;
            end
            [~,ind_lambda_opt] = min(mse);
            lambda_opt = val_lambda(ind_lambda_opt);
    end

end

