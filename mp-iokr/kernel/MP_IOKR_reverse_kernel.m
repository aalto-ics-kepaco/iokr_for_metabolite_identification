function [ score ] = MP_IOKR_reverse_kernel( KX_list, train_set, test_set, Y_train, Y_C, val_gamma, val_lambda, param )
%======================================================
% DESCRIPTION:
% MP-IOKR with reverse IOKR in the case of a kernel represention in output
%
% INPUTS:
% KX_list:      list of input kernel matrices (each of size n*n)
% train_set:    vector containing the training indices (length n_train)
% test_test:    vector containing the test indices (length n_test)
% Y_train:      matrix of size d*n_train containing the training output feature vectors
% Y_C:          Candidate sets
% val_gamma:    vector of strictly positive values among which the regularization 
%               parameter of reverse IOKR will be selected
% val_lambda:   vector of strictly positive values among which the regularization 
%               parameter of MP-IOKR will be selected
% param:        structure containing the MP-IOKR parameters
%   param.center:       binary value indicating if the input and output
%                       kernel/feature vectors should be centered (1) or not (0)
%   param.mkl:          string indicating the MKL algorithm used for kernel combination 
%                       ('alignf' or 'unimkl')
%   param.rev_iokr:     string indicating if the input feature maps of the different
%                       input kernels should be learnt jointly ('joint') or separately ('separate')
%
% OUTPUTS:
% score:        cell of length n_test, each entry being a vector containing 
%               the scores associated with the candidates in the corresponding
%               candidate set
%
%======================================================

    ker_center = param.center; % centering option

    % Learning kernel weights with Multiple Kernel Learning  
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
            for k = 1:n_k
                [KX_train_list{k}, ~] = input_kernel_center_norm(KX_list{k}, train_set, test_set, mp_iokr_param.center);
                KX_train_list{k} = w(k) * KX_train_list{k};
            end
    end
    
    % Build output kernel and compute a struct object containing the
    % kernel parameters and the informations needed for kernel centering and normalization
    KY_train = build_kernel(Y_train, Y_train, KY_par);
    OK.mean = mean(KY_train,1);
    
    KY_train_c = center(KY_train, OK.mean, ker_center); % centering
    OK.diag = diag(KY_train_c);
    
    OK.param = KY_par;
    OK.center = param.center;
    
    % Parameter selection
    % to do
    
    % Training the reverse IOKR model
    M = Train_reverse_IOKR_kernel(normmat(KY_train_c), gamma_opt);
    
    % Training the MP-IOKR model
    [Mean_KY_train_C, Cov_KY_train_C] = compute_cov_mean_kernel( Y_train, Y_C(train_set), OK);
    
    C = Train_MP_IOKR_reverse_kernel(KX_train, M, Mean_KY_train_C, Cov_KY_train_C, lambda);
    B = C \ KX_train_test;
    
    % Pre-image                                                  
    score = Preimage_MP_IOKR_reverse_kernel(M, Y_train, OK, Y_C(train_set), Y_C(test_set), B, Mean_KY_train_C);

end

