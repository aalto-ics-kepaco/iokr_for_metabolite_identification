function [ mse ] = MPIOKR_kernel_eval_mse( KX_list, Y, Y_C, KY_par, gamma_opt, ...
    mp_iokr_param, opt_param, debug_param )
%======================================================
% DESCRIPTION:
% Computation of the mean squared errors (mse) for different values of the 
% regularization parameter( kernel represention in output)
%
% INPUTS:
% KX_list:      cell array containing the training input kernel matrices
% Y:            matrix containing the training output vectors
% Y_C:          training candidates encapsulated in a CandidateSet object
% KY_par:       output kernel parameters
% mp_iokr_param:1*1 struct array containing information relative to
%               centering and multiple kernel learning
% select_param: 1*1 struct array containing information related to the parameter selection 
%
% OUTPUT:
% mse:          vector containing the mse obtained for each regularization parameter
%               in select_param.lambda
%
%====================================================== 

    % Vector containing the possible values for the regularization parameters
    val_lambda = opt_param.val_lambda;
    
    n_kx = length(KX_list);
    n = size(Y,2);
    
    [KY_tilde, KY_S_C, V, D] = build_tilde_kernel (Y, Y_C, KY_par, ...
        mp_iokr_param, debug_param);
    KY_S = output_kernel_preprocessing_train (Y, KY_par, ...
        mp_iokr_param.center);
    
    % Parameter selection using inner cross-validation
    c = opt_param.cv_partition; % CV partition
    n_folds = c.NumTestSets; % number of folds
    
    mse_cv = zeros(length(val_lambda), n_folds);
    
    for j = 1:n_folds % Cross-validation
        
        if (debug_param.verbose)
            tic;
            fprintf ('Inner fold: %d/%d\n', j, n_folds);
        end % if
        
        ind_train_S = find(training(c,j))'; % indices of the training examples
        ind_test_S = find(test(c,j))'; % indices of the test examples

        n_train = length(ind_train_S);
        n_test = length(ind_test_S);
       
        ind_train_C = find(sum(V(:,ind_train_S),2))'; % indices of the candidate sets of the training examples
        ind_test_C = find(sum(V(:,ind_test_S),2))'; % indices of the candidate sets of the test examples
        
        n_train_C = length(ind_train_C);
        n_test_C = length(ind_test_C);
        
        ind_train_tilde = [ind_train_S,ind_train_C+n]; % indices corresponding to the training and training candidates examples
        ind_test_tilde = [ind_test_S,ind_test_C+n]; % indices corresponding to the test and test candidates examples
        
        % input kernels
        KX_list_train = cellfun(@(x) x(ind_train_S,ind_train_S), KX_list, 'UniformOutput', false);
        KX_list_train_test = cellfun(@(x) x(ind_train_S,ind_test_S), KX_list, 'UniformOutput', false);

        KX_train = blkdiag(KX_list_train{1:n_kx});
        KX_train_test = cell2mat(KX_list_train_test);
        clear KX_train_test_list
        
        % output kernels
        KY_tilde_train = KY_tilde(ind_train_tilde,ind_train_tilde);
        KY_tilde_test = KY_tilde(ind_test_tilde,ind_test_tilde);
        KY_tilde_train_test = KY_tilde(ind_train_tilde,ind_test_tilde);
        
        KY_S_C_train = KY_S_C(ind_train_S,ind_train_C);
        KY_S_C_train_test = KY_S_C(ind_train_S,ind_test_C);
        
        V_train = V(ind_train_C,ind_train_S);
        V_test = V(ind_test_C,ind_test_S);
        
        D_train = D(ind_train_C,ind_train_C);
        D_test = D(ind_test_C,ind_test_C);
        
        KY_S_train = KY_S(ind_train_S, ind_train_S);
        M_train_test_c = cell(n_kx,1);
        M_train_c      = cell(n_kx,1);
        for k = 1:n_kx
            M_train_test_c{k} = (gamma_opt(k) * eye(n_train) + KY_S_train) \ KY_S_C_train_test;
            M_train_c{k}      = (gamma_opt(k) * eye(n_train) + KY_S_train) \ KY_S_C_train;
        end
        M_train_test = cell2mat (M_train_test_c);
        M_train      = cell2mat (M_train_c);
        clear M_train_test_c M_train_c;
        
%         KX_train_test_tilde = [KX_train_test - KX_train*M_train_test*D_test.^2*V_test, ...
%                                KX_train*M_train_test*(eye(n_test_C) - D_test.^2*(V_test*V_test'))*D_test];
        KX_train_test_tilde = [KX_train_test - KX_train*M_train_test*D_test^2*V_test, ...
                               KX_train*M_train_test*(eye(n_test_C) - D_test^2*(V_test*V_test'))*D_test];
                           
%         fprintf ('Mean diag KY_tilde: %f\n', mean (diag (KY_tilde)));
        % mean (diag (KX_train_test_tilde)); 
                                   
        Ic = cell(n_kx,1);
        for k = 1:n_kx
            Ic{k} = eye(n_train);
        end
        I = cell2mat(Ic);
        clear Ic;

        % Training MP-IOKR with the reverse IOKR
%         A = [I - M_train*D_train.^2*V_train, M_train*(eye(n_train_C) - D_train.^2*(V_train*V_train'))*D_train];
        A = [I - M_train*D_train^2*V_train, M_train*(eye(n_train_C) - D_train^2*(V_train*V_train'))*D_train];
        AAt = A*A';

        for il = 1:length(val_lambda)
            lambda = val_lambda(il);
            
%             if (debug_param.verbose)
%                 fprintf ('Lambda: %e\n', lambda);
%             end % if

            % Training
            C = (lambda * eye(n_kx*n_train) + KX_train * AAt);
            B = C \ KX_train_test_tilde;
            
            mse_cv(il,j) = 1/n_test*trace(B'*A*KY_tilde_train*A'*B + KY_tilde_test - 2*KY_tilde_train_test'*A'*B);
        end
        
        if (debug_param.verbose)
            toc;
        end % if
    end
    mse = mean(mse_cv, 2); % mean over the different folds

end

