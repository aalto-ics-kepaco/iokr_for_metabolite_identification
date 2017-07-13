function [ mse ] = MPIOKR_feature_eval_mse( KX_list, Psi, Y_C, M, mp_iokr_param, select_param )
%======================================================
% DESCRIPTION:
% Computation of the mean squared errors (mse) for different values of the 
% regularization parameter
%
% INPUTS:
% KX_list:          cell array containing the training input kernel matrices
% Psi:              training output feature vectors
% Y_C:              cell array containing the training candidate output vectors
% M
% mp_iokr_param:    1*1 struct array containing information relative to
%                   centering and multiple kernel learning
% select_param:     1*1 struct array containing information related to the parameter selection
%
% OUTPUT:
% mse:          vector containing the mse obtained for each regularization parameter
%               in select_param.lambda
%
%====================================================== 

    % Vector containing the possible values for the regularization parameter
    val_lambda = select_param.lambda;
    
    n_kx = length(KX_list);
    
    [d,~] = size(Psi);

    % MP-IOKR: parameter selection using inner CV

    c = select_param.cv_partition; % CV partition
    n_folds = select_param.num_folds; % number of folds

    mse_cv = zeros(length(val_lambda), n_folds);

    for j = 1:n_folds % Cross-validation
        
        train_set = find(training(c,j));
        test_set = find(test(c,j));
        
        n_train = length(train_set);
        n_test = length(test_set);

        KX_list_train = cellfun(@(x) x(train_set,train_set), KX_list, 'UniformOutput', false);
        KX_list_train_test = cellfun(@(x) x(train_set,test_set), KX_list, 'UniformOutput', false);

        KX_train = blkdiag(KX_list_train{:});
        KX_train_test = cell2mat(KX_list_train_test);
        clear KX_train_test_cv_list

        Y_C_train = Y_C.getSubset(train_set);
        Y_C_test = Y_C.getSubset(test_set);

        [Mean_Psi_C_train, Cov_Psi_C_train] = Compute_cov_mean_feat(Y_C_train, mean_Y_train, mp_iokr_param.center, debug_param.verbose);
        [Mean_Psi_C_test, Cov_Psi_C_test] = Compute_cov_mean_feat(Y_C_test, mean_Y_train, mp_iokr_param.center, debug_param.verbose);
        
        Ic = cell(n_kx,1);
        for k = 1:n_kx
            Ic{k} = eye(n_train);
        end
        I = cell2mat(Ic);
        clear Ic;
        
        A1 = I - M * Mean_Psi_C_train;
        PsiAt = (Psi(:,train_set) - Mean_Psi_C_train) * A1' + Cov_Psi_C_train*M';
        AAt = A1 * A1' + M * Cov_Psi_C_train * M';


        for il = 1:length(val_lambda)
            lambda = val_lambda(il);

            C = PsiAt / (lambda * eye(n_kx*n_train) + KX_train * AAt);

            % Computation of the MP-error
            P = C * KX_train * M;

            mse_cv(il,j) = 1/n_test*norm((C*KX_train_test - P*Mean_Psi_C_test) - (Psi(:,test_set) - Mean_Psi_C_test),'fro')...
                + trace((P - eye(d))*Cov_Psi_C_test*(P - eye(d))');

        end
    end
    mse = mean(mse_cv, 2); % mean over the different folds
            


end

