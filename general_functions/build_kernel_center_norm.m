function [ K_nc ] = build_kernel_center_norm( X1, X2, mean_K1, D1_c, kernel_opt, ker_center )
%======================================================
% DESCRIPTION:
% Compute the gram matrix between the sets X1 and X2 and add centering and normalization
%
% INPUTS:
% X1:           matrix of size d * n1
% X2:           matrix of size d * n2
% mean_K1:      mean of the Gram matrix of the kernel between the set X1 and itself (only use if X2 ~= X1)
% D1_c          diagnoal vector of K1_c
% kernel_opt:   structure containing the kernel type and parameter(s). It can
%               be defined using the set_kernel_opt.m function
% ker_center:   binary value indicating if the Gram matrix has to be centered or not
%
% OUTPUTS:
% K_nc:           Gram matrix of size n1 * n2
%
%======================================================
    
    K = build_kernel(X1, X2, kernel_opt);
    K2 = build_kernel(X2, X2, kernel_opt);
    
    mean_K = mean(K,1);

    % Centering
    if ker_center == 1
        K_c = center(K, mean_K1, ker_center, mean_K1', mean_K);
        K2_c = center(K2, mean_K1, ker_center, mean_K', mean_K);
    else
        K_c = K;
        K2_c = K2;
    end
    
    % Normalization
    K_nc = normmat(K_c, D1_c, diag(K2_c));

end

