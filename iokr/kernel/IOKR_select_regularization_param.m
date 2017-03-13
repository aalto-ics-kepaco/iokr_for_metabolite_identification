function [ mse ] = IOKR_select_regularization_param( KX_train, KY_train, select_param )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    val_lambda = select_param.lambda;

    switch select_param.cv_type
        
        % Parameter selection with inner cross-validation
        case 'cv'
        
            c = select_param.cv_partition;
            n_folds = select_param.num_folds;

            mse_cv = zeros(1, length(val_lambda), n_folds);

            for j = 1:n_folds
                train_set_cv = find(training(c,j));
                test_set_cv = find(test(c,j));
                
                n_train_cv = length(train_set_cv);
                n_test_cv = length(test_set_cv);

                KX_train_cv = KX_train(train_set_cv, train_set_cv);
                KX_train_test_cv = KX_train(train_set_cv, test_set_cv);

                KY_train_cv = KY_train(train_set_cv, train_set_cv);
                KY_train_test_cv = KY_train(train_set_cv, test_set_cv);

                for il = 1:length(val_lambda)

                    % Training and prediction                    
                    C_cv = val_lambda(il)*eye(n_train_cv) + KX_train_cv;
               
                    B_cv = C_cv \ KX_train_test_cv;

                    % Compute the mean squared error
                    mse_cv(1,il,j) = 1 + 1/n_test_cv*...
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
                C = val_lambda(il)*eye(n_train) + KX_train;

                % Prediction on the training set
                B = C \ KX_train;

                % Compute the mean squared error
                LOOE = (eye(n_train)-B) / diag(diag(eye(n_train)-B));
                mse(1,il) = 1/n_train * trace(LOOE' * KY_train * LOOE);
            end
    end

end

