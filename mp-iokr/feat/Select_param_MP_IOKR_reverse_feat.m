function [ lambda_opt, mp_err ] = Select_param_MP_IOKR_reverse_feat( ...
    KX_train_list, Y_train, Y_C_train, gamma_opt, ...
    opt_param, mp_iokr_param, data_param, debug_param)
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
    if (debug_param.verbose)
        % Create some stop-watches
        sw_Train_reverse_IOKR_feat    = StopWatch ('Train_reverse_IOKR_feat (g-function)');
        sw_compute_cov_mean_feat      = StopWatch ('Compute_cov_mean_feat (train & test)');
        sw_train_MP_IOKR_reverse_feat = StopWatch ('Train_MP_IOKR_reverse_feat (h-function)');
        sw_calculate_mp_error         = StopWatch ('Calculate mp-error');
        sw_find_hyper_parameter       = StopWatch ('Find hyper-parameter (Train_reverse_IOKR_feat_reverse_feat & Calculate mp-error)');
        sw_whole_loop                 = StopWatch ('WHOLE LOOP');
    end % if 

    [d,n_train] = size(Y_train);
      
    % Training output feature vectors
    mean_Y_train = mean(Y_train,2);
    Psi_train = norma(Y_train, mean_Y_train, mp_iokr_param.center);
    
    % Selection of the regularization parameter of MP-IOKR using an inner
    % cross-validation experiment
    if (data_param.usePreCalcStat)
        n_folds = data_param.cv.NumTestSets;
        cv = data_param.cv;
    else
        n_folds = opt_param.nInnerFolds;
        cv = cvpartition(n_train, 'k', n_folds);
    end % if
            
    mp_err = zeros(length(opt_param.val_lambda), n_folds);
    for foldIdx = 1:n_folds     
        if (debug_param.verbose) ; sw_whole_loop.start() ; end % if
            
        fprintf ('Inner fold: %d/%d\n', foldIdx, n_folds);
        
        % Defining training and test sets       
        % Insert some assertions to prevent problems with pre-calulated
        % data. Can be removed later.
        assert (numel (training_my(cv,foldIdx)) == size (KX_train_list{1}, 1), ...
            'Upps?!: Lenght of the binary training selection vector should match the dimension of the kernel. usePreCalcStat = %d', ...
            data_param.usePreCalcStat);
        assert (numel (test_my(cv,foldIdx)) == size (KX_train_list{1}, 1), ...
            'Upps?!: Lenght of the binary test selection vector should match the dimension of the kernel. usePreCalcStat = %d', ...
            data_param.usePreCalcStat);
        
        train_set_cv = find(training_my(cv,foldIdx));        
        test_set_cv = find(test_my(cv,foldIdx));
        
        KX_train_cv_list = cellfun(@(x) x(train_set_cv,train_set_cv), KX_train_list, 'UniformOutput', false);
        KX_train_test_cv_list = cellfun(@(x) x(train_set_cv,test_set_cv), KX_train_list, 'UniformOutput', false);
        
        KX_train_test_cv = cell2mat(KX_train_test_cv_list);
        clear KX_train_test_cv_list
        
        Psi_train_cv = Psi_train(:,train_set_cv);
        Psi_test_cv = Psi_train(:,test_set_cv);
        
        % Train the reverse IOKR model
        if (debug_param.verbose) ; sw_Train_reverse_IOKR_feat.start() ; end % if      
        
        M_cv = Train_reverse_IOKR_feat(Psi_train_cv, gamma_opt);        
        
        if (debug_param.verbose) ; sw_Train_reverse_IOKR_feat.stop() ; end % if
             
        % Compute the mean and covariance of the candidate output feature
        % vectors
        if (data_param.usePreCalcStat)             
            % Train
            Mean_Psi_C_train_cv = data_param.stats_cv(foldIdx, 1).Mean_Psi_C_train_cv;
            Cov_Psi_C_train_cv = data_param.stats_cv(foldIdx, 1).Cov_Psi_C_train_cv;
            % Test
            Mean_Psi_C_test_cv = data_param.stats_cv(foldIdx, 1).Mean_Psi_C_test_cv;
            Cov_Psi_C_test_cv = data_param.stats_cv(foldIdx, 1).Cov_Psi_C_test_cv;
            
            clear stats_cv;
        else
            assert (Y_C_train.getNumberOfExamples() == numel (training_my(cv,foldIdx)), 'Upps?!')
            assert (Y_C_train.getNumberOfExamples() == numel (test_my(cv,foldIdx)), 'Upps?!')
            
            Y_C_train_cv = Y_C_train.getSubset (train_set_cv);
            Y_C_test_cv = Y_C_train.getSubset (test_set_cv);
        
            if (debug_param.verbose) ; sw_compute_cov_mean_feat.start(); end % if
            
            [Mean_Psi_C_train_cv, Cov_Psi_C_train_cv] = Compute_cov_mean_feat(Y_C_train_cv, mean_Y_train, mp_iokr_param.center, debug_param.verbose);
            [Mean_Psi_C_test_cv, Cov_Psi_C_test_cv] = Compute_cov_mean_feat(Y_C_test_cv, mean_Y_train, mp_iokr_param.center, debug_param.verbose);
            
            if (debug_param.verbose) ; sw_compute_cov_mean_feat.stop() ; end % if
        end % if
        
        if (debug_param.verbose) ; sw_find_hyper_parameter.start(); end % if
        
        for val_lambda_idx = 1:length(opt_param.val_lambda)
            lambda = opt_param.val_lambda(val_lambda_idx);
            
            % Training MP-IOKR
            if (debug_param.verbose) ; sw_train_MP_IOKR_reverse_feat.start(); end % if
            
            C_cv = Train_MP_IOKR_reverse_feat(KX_train_cv_list, Psi_train_cv, M_cv, Mean_Psi_C_train_cv, Cov_Psi_C_train_cv, lambda);
            
            if (debug_param.verbose) ; sw_train_MP_IOKR_reverse_feat.stop() ; end % if
            
            % Computation of the MP-error
            if (debug_param.verbose) ; sw_calculate_mp_error.start(); end % if

            KX_train_cv = blkdiag (KX_train_cv_list{:});
            P = C_cv*KX_train_cv*M_cv;
            mp_err(val_lambda_idx,foldIdx) = norm((C_cv*KX_train_test_cv - P*Mean_Psi_C_test_cv) - (Psi_test_cv - Mean_Psi_C_test_cv),'fro')...
                + trace((P - eye(d))*Cov_Psi_C_test_cv*(P - eye(d))');
            mp_err(val_lambda_idx,foldIdx) = 1/length(test_set_cv) * mp_err(val_lambda_idx,foldIdx);
            
            if (debug_param.verbose) ; sw_calculate_mp_error.stop() ; end % if
        end
        
        if (debug_param.verbose) ; sw_find_hyper_parameter.stop(); end % if
        
        if (debug_param.verbose)
            sw_whole_loop.stop();
                        
            sw_Train_reverse_IOKR_feat.showAvgTime();
            sw_compute_cov_mean_feat.showAvgTime();
            sw_find_hyper_parameter.showAvgTime();
            fprintf ('  > '); sw_train_MP_IOKR_reverse_feat.showAvgTime();
            fprintf ('  > '); sw_calculate_mp_error.showAvgTime();
            
            fprintf ('----------\n')
            
            avgTime_sumOfWatches = sw_Train_reverse_IOKR_feat.avgTime_ + ...
                sw_compute_cov_mean_feat.avgTime_ + sw_find_hyper_parameter.avgTime_;
            fprintf ('%.3fs --- ALL stop-watches\n', avgTime_sumOfWatches);
            sw_whole_loop.showAvgTime();
            
            fprintf ('----------\n')
            
            fprintf ('%.3fs --- OVERHEAD\n', sw_whole_loop.avgTime_ - avgTime_sumOfWatches);
        end % if
    end
    
    [~, ind_lambda_opt] = min(mean(mp_err,2));
    lambda_opt = opt_param.val_lambda(ind_lambda_opt);
end

