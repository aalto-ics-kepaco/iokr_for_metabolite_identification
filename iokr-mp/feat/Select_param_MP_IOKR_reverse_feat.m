function [ gamma_opt, lambda_opt ] = Select_param_MP_IOKR_reverse_feat( ...
    KX_train_list, Y_train, Y_C_train, opt_param, mp_iokr_param, data_param)
%======================================================
% DESCRIPTION:
% Hyperparameter selection for MP-IOKR with reverse IOKR in the case of a feature
% represention in output
%
% INPUTS:
% KX_train_list:    list of input kernel matrices on the training set (size n_train*n_train)
% Y_train:          matrix containing the training output vectors (size d*n_train)
% Y_C_train:        training candidate sets
% opt_param:    structure containing the optimization parameters
%   opt_param.val_gamma:  vector of strictly positive values among which 
%                         the regularization parameter of reverse IOKR will 
%                         be selected
%   opt_param.val_lambda: vector of strictly positive values among which 
%                         the regularization  parameter of MP-IOKR will be 
%                         selected
% mp_iokr_param:    structure containing the MP-IOKR parameters
%   mp_iokr_param.center:   binary value indicating if the input and output
%                           kernel/feature vectors should be centered (1) or not (0)
%   mp_iokr_param.mkl:      string indicating the MKL algorithm used for 
%                           kernel combination ('alignf' or 'unimkl')
%   mp_iokr_param.rev_iokr: string indicating if the input feature maps of 
%                           the different input kernels should be learnt 
%                           jointly ('joint') or separately ('separate')
% data_param:   structure containing the parameters for the data handling
%   data_param.usePreCalcStat:  binary indicating whether pre-calculated
%                               candidate-set mean vectors and covariance
%                               matrices should be used
%   ONLY needed (and accessed) if data_param.usePreCalcStat is to true.
%   data_param.cv               pre-fixed cross-validation folds for the
%                               parameter selection 
%   data_param.statsMatObj      Object of class matlab.io.MatFile storing
%                               the statistics of the candidate sets for
%                               corresponding to each fold. 
% OUTPUTS:
% gamma_opt, lambda_opt: selected regularization parameters
%
%======================================================
    [d,n_train] = size(Y_train);
      
    % Training output feature vectors
    mean_Y_train = mean(Y_train,2);
    Psi_train = norma(Y_train, mean_Y_train, mp_iokr_param.center);
    
    % Selection of the regularization parameter(s) of reverse IOKR
    gamma_opt = Select_param_reverse_IOKR(KX_train_list, Psi_train, opt_param.val_gamma);
    
    % Selection of the regularization parameter of MP-IOKR using an inner
    % cross-validation experiment
    if (data_param.usePreCalcData)
        n_folds = data_param.cv.NumTestTests;
        c = opt_param.cv;
    else
        n_folds = 10;
        c = cvpartition(n_train, 'k', n_folds);
    end % if
    
    e = zeros(length(opt_param.val_lambda), n_folds);
    for i = 1:n_folds       
        % Defining training and test sets       
        % Insert some assertions to prevent problems with pre-calulated
        % data. Can be removed later.
        assert (numel (training_my(c,i)) == size (KX_train_list{1}, 1), ...
            'Upps?!: Lenght of the binary training selection vector should match the dimension of the kernel. usePreCalcStat = %d', ...
            data_param.usePreCalcData);
        assert (numel (test_my(c,i)) == size (KX_train_list{1}, 1), ...
            'Upps?!: Lenght of the binary test selection vector should match the dimension of the kernel. usePreCalcStat = %d', ...
            data_param.usePreCalcData);
        
        train_set_cv = find(training_my(c,i));        
        test_set_cv = find(test_my(c,i));
        
        KX_train_cv_list = cellfun(@(x) x(train_set_cv,train_set_cv), KX_train_list, 'UniformOutput', false);
        KX_train_test_cv_list = cellfun(@(x) x(train_set_cv,test_set_cv), KX_train_list, 'UniformOutput', false);
        
        KX_train_test_cv = cell2mat(KX_train_test_cv_list);
        clear KX_train_test_cv_list
        
        Psi_train_cv = Psi_train(:,train_set_cv);
        Psi_test_cv = Psi_train(:,test_set_cv);
        
        % Train the reverse IOKR model
        M_cv = Train_reverse_IOKR(Psi_train_cv, gamma_opt);
             
        % Compute the mean and covariance of the candidate output feature
        % vectors
        if (data_param.usePreCalcData) 
            stats_cv = data_param.statMatObj.stats_cv(i, 1);
            
            % Train
            Mean_Psi_C_train_cv = stats_cv.Mean_Psi_C_train_cv;
            Cov_Psi_C_train_cv = stast.Cov_Psi_C_train_cv;
            % Test
            Mean_Psi_C_test_cv = stats_cv.Mean_Psi_C_test_cv;
            Cov_Psi_C_test_cv = stast.Cov_Psi_C_test_cv;
            
            clear stats_cv;
        else
            assert (Y_C_train.getNumberOfExamples() == numel (train_set_cv_logical), 'Upps?!')
            assert (Y_C_train.getNumberOfExamples() == numel (test_set_cv_logical), 'Upps?!')
            
            Y_C_train_cv = Y_C_train.getSubset (train_set_cv);
            Y_C_test_cv = Y_C_train.getSubset (test_set_cv);

            [Mean_Psi_C_train_cv, Cov_Psi_C_train_cv] = Compute_cov_mean_feat(Y_C_train_cv, mean_Y_train, mp_iokr_param.center);
            [Mean_Psi_C_test_cv, Cov_Psi_C_test_cv] = Compute_cov_mean_feat(Y_C_test_cv, mean_Y_train, mp_iokr_param.center);
        end % if
        
        for j = 1:length(opt_param.val_lambda)
            lambda = opt_param.val_lambda(j);
            
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
    lambda_opt = opt_param.val_lambda(ind_lambda_opt);
end

