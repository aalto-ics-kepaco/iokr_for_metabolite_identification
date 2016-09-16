function [ gamma_opt, lambda_opt ] = Select_param_MP_IOKR_reverse_feat( KX_train_list, Y_train, Y_C_train, val_gamma, val_lambda, ker_center )
%======================================================
% DESCRIPTION:
% Hyperparameter selection for MP-IOKR with reverse IOKR in the case of a feature
% represention in output
%
% INPUTS:
% KX_train_list:    list of input kernel matrices on the training set (size n_train*n_train)
% Y_train:          matrix containing the training output vectors (size d*n_train)
% Y_C_train:        training candidate sets
% val_gamma:        vector of strictly positive values among which the regularization
%                   parameter of reverse IOKR will be selected
% val_lambda:       vector of strictly positive values among which the regularization
%                   parameter of MP-IOKR will be selected
% ker_center:       value of 0 of 1 indicating if the output feature
%                   vectors should be centered.
%
% OUTPUTS:
% gamma_opt, lambda_opt: selected regularization parameters
%
%======================================================

    [d,n_train] = size(Y_train);
    
    n_kx = length(KX_train_list); % number of input kernels
    
    % Training output feature vectors
    mean_Y_train = mean(Y_train,2);
    Psi_train = norma(Y_train, mean_Y_train, ker_center);
    
    % Selection of the regularization parameter(s) of reverse IOKR
    gamma_opt = Select_param_reverse_IOKR(KX_train_list, Psi_train, val_gamma);
    
    % Selection of the regularization parameter of MP-IOKR using an inner
    % cross-validation experiment
    n_folds = 10;
    c = cvpartition(n_train, 'k', n_folds);
    e = zeros(length(val_lambda), n_folds);
    for i = 1:n_folds
        
        % Defining training and test sets
        train_set_cv = find(training(c,i));
        test_set_cv = find(test(c,i));
        
        KX_train_cv_list = cellfun(@(x) x(train_set_cv,train_set_cv), KX_train_list, 'UniformOutput', false);
        KX_train_test_cv_list = cellfun(@(x) x(train_set_cv,test_set_cv), KX_train_list, 'UniformOutput', false);
        
        KX_train_test_cv = cell2mat(KX_train_test_cv_list);
        clear KX_train_test_cv_list
        
        Psi_train_cv = Psi_train(:,train_set_cv);
        Psi_test_cv = Psi_train(:,test_set_cv);
        
        Y_C_train_cv = Y_C_train(train_set_cv);
        Y_C_test_cv = Y_C_train(test_set_cv);
        
        % Train the reverse IOKR model
        M_cv = Train_reverse_IOKR(Psi_train_cv, gamma_opt);
             
        % Compute the mean and covariance of the candidate output feature
        % vectors
        [Mean_Psi_C_train_cv, Cov_Psi_C_train_cv] = Compute_cov_mean_feat(Y_C_train_cv, mean_Y_train, ker_center);
        [Mean_Psi_C_test_cv, Cov_Psi_C_test_cv] = Compute_cov_mean_feat(Y_C_test_cv, mean_Y_train, ker_center);
        
        for j = 1:length(val_lambda)
            lambda = val_lambda(j);
            
            % Training MP-IOKR
            C_cv = Train_MP_IOKR_reverse_feat(KX_train_cv_list, Psi_train_cv, M_cv, Mean_Psi_C_train_cv, Cov_Psi_C_train_cv, lambda);
            
            % Computation of the MP-error
            P = C_cv*KX_train_cv*M_cv;
            e(j,i) = norm((C_cv*KX_train_test_cv - P*Mean_Psi_C_test_cv) - (Psi_test_cv - Mean_Psi_C_test_cv),'fro')...
                + trace((P - eye(d))*Cov_Psi_C_test_cv*(P - eye(d))');
            e(j,i) = 1/length(test_set_cv) * e(j,i);
        end
    end
    [~, ind_lambda_opt] = min(mean(e,2));
    lambda_opt = val_lambda(ind_lambda_opt);

end

