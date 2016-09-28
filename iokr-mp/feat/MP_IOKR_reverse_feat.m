function [ score, debug_info ] = MP_IOKR_reverse_feat(KX_list, Y_train, Y_C, ...
    opt_param, mp_iokr_param, data_param, debug_param)
%======================================================
% DESCRIPTION:
% MP-IOKR with reverse IOKR in the case of a feature represention in output
%
% INPUTS:
% KX_list:      list of input kernel matrices (each of size n*n)
% Y_train:      matrix of size d*n_train containing the training output feature vectors
% Y_C:          Candidate sets
% opt_param:    structure containing the optimization parameters
%
%   opt_param.val_gamma:  vector of strictly positive values among which 
%                         the regularization parameter of reverse IOKR will 
%                         be selected
%   opt_param.val_lambda: vector of strictly positive values among which 
%                         the regularization  parameter of MP-IOKR will be 
%                         selected
%
% mp_iokr_param:    structure containing the MP-IOKR parameters
%
%   mp_iokr_param.center:   binary value indicating if the input and output
%                           kernel/feature vectors should be centered (1) or not (0)
%   mp_iokr_param.mkl:      string indicating the MKL algorithm used for 
%                           kernel combination ('alignf' or 'unimkl')
%   mp_iokr_param.rev_iokr: string indicating if the input feature maps of 
%                           the different input kernels should be learnt 
%                           jointly ('joint') or separately ('separate')
%
% data_param:   structure containing the parameters for the data handling
%
%   data_param.train_set        binary vector containing indicating the
%                               examples used for training (length n, 
%                               sum (data_param.train_set) == n_train)
%   data_param.test_set         binary vector containing indicating the
%                               examples used for testing (length n, 
%                               sum (data_param.test_set) == n_test)
%   data_param.usePreCalcStat:  binary indicating whether pre-calculated
%                               candidate-set mean vectors and covariance
%                               matrices should be used
%   ONLY needed (and accessed) if DATA_PARAM.USEPRECALCSTAT is to true.
%   data_param.cv               pre-fixed cross-validation folds for the
%                               parameter selection 
%                               This overwrites DATA_PARAM.TRAIN_SET and 
%                               DATA_PARAM.TEST_SET
%   data_param.stats      Object of class matlab.io.MatFile storing
%                               the statistics of the candidate sets for
%                               corresponding to each fold. 
%
% OUTPUTS:
% score:        cell of length n_test, each entry being a vector containing 
%               the scores associated with the candidates in the corresponding
%               candidate set
%
%======================================================
    if (debug_param.verbose)
        sw_mkl_weights                     = StopWatch ('mkl-weights');
        sw_input_kernel_processing         = StopWatch ('input-kernel-processing');
        sw_train_reverse_IOKR              = StopWatch ('Train_reverse_IOKR');
        sw_compute_cov_mean_feat           = StopWatch ('Compute_cov_mean_feat');
        sw_train_MP_IOKR_reverse_feat      = StopWatch ('Train_MP_IOKR_reverse_feat');
        sw_prediction_MP_IOKR_reverse_feat = StopWatch ('Prediction_MP_IOKR_reverse_feat');
        sw_preimage_MP_IOKR_feat           = StopWatch ('Preimage_MP_IOKR_feat');
    end % if

    train_set = find (data_param.train_set);
    test_set  = find (data_param.test_set);

    KX_list_train = cellfun(@(x) x(train_set,train_set), KX_list, 'UniformOutput', false);
    
    %% Learning kernel weights with Multiple Kernel Learning  
    if (debug_param.verbose) ; sw_mkl_weights.start() ; end % if
    
    w = mkl_weight(mp_iokr_param.mkl, KX_list_train, normmat(Y_train'*Y_train));
    
    if (debug_param.verbose)
        sw_mkl_weights.stop();
        sw_mkl_weights.showAvgTime();
    end % if
        
    clear KX_list_train
    
    %% Centering and normalization of input kernel
    if (debug_param.verbose) ; sw_input_kernel_processing.start() ; end % if
    
    n_kx = length(KX_list);
    switch mp_iokr_param.rev_iokr
        case 'joint' 
            % Computation of the combined input kernel
            [KX_train_combined, KX_train_test_combined] = mkl_combine_train_test(KX_list, train_set, test_set, w, mp_iokr_param.center);
            KX_train_list = {KX_train_combined};
            KX_train_test_list = {KX_train_test_combined};
            clear KX_train_combined KX_train_test_combined
        case 'separate'
            KX_train_list = cell(n_kx,1);
            KX_train_test_list = cell(n_kx,1);
            for k = 1:n_kx
                [KX_train_list{k}, KX_train_test_list{k}] = input_kernel_center_norm(KX_list{k}, train_set, test_set, mp_iokr_param.center);
                KX_train_list{k} = w(k) * KX_train_list{k};
                KX_train_test_list{k} = w(k) * KX_train_test_list{k};
            end
    end
    
    if (debug_param.verbose)
        sw_input_kernel_processing.stop();
        sw_input_kernel_processing.showAvgTime();
    end % if

    %% Training output feature vectors
    mean_Y_train = mean(Y_train,2);
    Psi_train = norma(Y_train, mean_Y_train, mp_iokr_param.center);
    
    if (data_param.usePreCalcStat)
        Y_C_train = []; % Y_C_train is never accessed in this case.
    else
        Y_C_train = Y_C.getSubset (train_set);
    end % if
    
    %% Parameter selection    
    fprintf ('Parameter selection\n');
    debug_info = struct ();
    
    [gamma_opt, lambda_opt, debug_info.mp_err] = Select_param_MP_IOKR_reverse_feat ( ...
        KX_train_list, Y_train, Y_C_train, opt_param, mp_iokr_param, data_param, debug_param);
%     lambda_opt = 0.5;
%     gamma_opt = 0.5;
%     debug_info.mp_err = [];
    
    debug_info.gamma_opt = gamma_opt;
    debug_info.lambda_opt = lambda_opt;
    
    fprintf ('Selected parameter: lambda = %e, gamma = %e\n', ...
        debug_info.lambda_opt, debug_info.gamma_opt);
    
    if (debug_param.verbose)
        fprintf ('MP-ERR: mean over folds (1 x n_folds) & (n_val_lambda x n_folds):\n');
        disp ([mean(debug_info.mp_err, 2), debug_info.mp_err]);
    end % if
    
    %% Training the reverse IOKR model
    if (debug_param.verbose) ; sw_train_reverse_IOKR.start() ; end % if
    
    M = Train_reverse_IOKR(Psi_train, gamma_opt);
    
    if (debug_param.verbose) 
        sw_train_reverse_IOKR.stop();
        sw_train_reverse_IOKR.showAvgTime();
    end % if
    
    %% Training the MP-IOKR model
    if (data_param.usePreCalcStat) 
        % Train
        Mean_Psi_C_train = data_param.Mean_Psi_C_train;
        Cov_Psi_C_train = data_param.Cov_Psi_C_train;
        
        clear stats;
    else
        if (debug_param.verbose) ; sw_compute_cov_mean_feat.start() ; end % if
            
        [Mean_Psi_C_train, Cov_Psi_C_train] = Compute_cov_mean_feat(Y_C_train, mean_Y_train, mp_iokr_param.center, debug_param.verbose);
        
        if (debug_param.verbose)
            sw_compute_cov_mean_feat.stop();
            sw_compute_cov_mean_feat.showAvgTime();
        end % if
    end % if
    if (debug_param.verbose) ; sw_train_MP_IOKR_reverse_feat.start() ; end % if
    
    C = Train_MP_IOKR_reverse_feat(KX_train_list, Psi_train, M, Mean_Psi_C_train, Cov_Psi_C_train, lambda_opt);
    
    if (debug_param.verbose)
        sw_train_MP_IOKR_reverse_feat.stop();
        sw_train_MP_IOKR_reverse_feat.showAvgTime();
    end % if
    %% Prediction on the test set  
    if (debug_param.verbose) ; sw_prediction_MP_IOKR_reverse_feat.start() ; end % if
    
    Psi_pred = Prediction_MP_IOKR_reverse_feat(KX_train_test_list, C);
    
    if (debug_param.verbose)
        sw_prediction_MP_IOKR_reverse_feat.stop();
        sw_prediction_MP_IOKR_reverse_feat.showAvgTime();
    end % if    
    %% Preimage
    Y_C_test = Y_C.getSubset (test_set);
    
    if (debug_param.verbose) ; sw_preimage_MP_IOKR_feat.start() ; end % if
    
    score = Preimage_MP_IOKR_feat(Psi_pred, Y_C_test, mean_Y_train, mp_iokr_param.center);
    
    if (debug_param.verbose)
        sw_preimage_MP_IOKR_feat.stop();
        sw_preimage_MP_IOKR_feat.showAvgTime();
    end % if    
end