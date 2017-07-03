function [ train_model ] = Train_MPIOKR( KX_list_train, Y_train, Y_C_train, ky_param, mp_iokr_param, select_param, debug_param )
%======================================================
% DESCRIPTION:
% Training step of Magnitude-preserving IOKR
%
% INPUTS:
% KX_list_train:    cell array containing the training input kernel matrices
% Y_train:          matrix of size d*n_train containing the training output vectors
% output_param:     1*1 struct array containing the information related to the outputs
% select_param:     1*1 struct array containing information related to the parameter selection 
% iokr_param:       1*1 struct array containing information relative to
%                   centering and multiple kernel learning
%
% OUTPUT:
% train_model:      training model
%
%======================================================  

    n_train = size(Y_C_train,2);

    % Selection of the regularization parameter and of the output kernel parameters(s)
    [lambda_opt, gamma_opt, KY_par_opt, w_opt] = Select_param_MPIOKR(KX_list_train, Y_train, Y_C_train, ky_param, select_param, mp_iokr_param);

    % Input kernel processing
    if (strcmp (mp_iokr_param.rev_iokr, 'separate'))
        [KX_train, process_input] = mpiokr_input_kernel_preprocessing_train(KX_list_train, w_opt, mp_iokr_param, gamma_opt);
        gamma_opt = unique(gamma_opt);
    else
        [KX_train, process_input] = mpiokr_input_kernel_preprocessing_train(KX_list_train, w_opt, mp_iokr_param);
    end
    n_kx = length(KX_train);
    
    % Output processing
    if strcmp(ky_param.type,'linear')
        [Psi_train, process_output] = output_feature_preprocessing_train(Y_train, mp_iokr_param.center);
    else
        [KY_train, process_output] = output_kernel_preprocessing_train(Y_train, ky_param_all_comb(ip), mp_iokr_param.center);
    end
    
    % Computation of the matrices M and I
    Mc = cell(n_kx,1);
    for k = 1:n_kx
        if strcmp(KY_par_opt.type,'linear')
            Mc{k} = (gamma_opt(k) * eye(n_train) + (Psi_train'*Psi_train)) \ (Psi_train');
        else
            Mc{k} = inv(gamma_opt(k) * eye(n_train) + KY_train);
        end
    end
    M = cell2mat(Mc);
    
    Ic = cell(n_kx,1);
    for k = 1:n_kx
        Ic{k} = eye(n_train);
    end
    I = cell2mat(Ic);
    clear Ic;
    
    % Training
    if strcmp(ky_param.type,'linear')
        [Mean_Psi_C_train, Cov_Psi_C_train] = Compute_cov_mean_feat(Y_C_train, process_output.mean, mp_iokr_param.center, debug_param.verbose);
        
        A1 = I - M * Mean_Psi_C_train;
        AAt = A1 * A1' + M * Cov_Psi_C_train * M';
    else
        [~, KY_S_C, V, D] = build_tilde_kernel(Y_train, Y_C_train, KY_par_opt, mp_iokr_param);
        A1 = I - M * KY_S_C * D^2 * V;
        A2 = M * KY_S_C * (eye(size(D)) - D^2*(V*V')) * D;
        AAt = A1 * A1' + A2 * A2';
    end
    
    C = (lambda_opt * eye(n_kx*n_train) + KX_train * AAt);
    
    train_model = struct('C',C,'process_input',process_input,'process_output',process_output,'KY_par',KY_par_opt);
    
end

