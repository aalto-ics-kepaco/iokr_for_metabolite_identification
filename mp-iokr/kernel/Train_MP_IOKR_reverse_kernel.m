function [ C ] = Train_MP_IOKR_reverse_kernel(KX_train_list, M, Mean_KY_train_C, Cov_KY_train_C, lambda)
%======================================================
% DESCRIPTION:
% Training of MP-IOKR with reverse IOKR in the case of a kernel
% represention in output
%
% INPUTS:
% KX_train_list:    list containing one combined input kernel matrix on the training set (joint reverse IOKR) 
%                   or n_kx input kernel matrices on the training set, each of them being multiplied by the 
%                   corresponding mkl weight (separate reverse IOKR for each kernel)
% M:                list of 1 (joint case)/n_kx (separate case) reverse IOKR model(s)
% Mean_KY_train_C:  matrix of size n_train*n_train containing the averaged 
%                   output kernel values between the training set and the candidates
% Cov_KY_train_C:   matrix of size n_train*n_train containing the covariance of the 
%                   output kernel values between the training set and the candidates
% lambda:           regularization parameter of MP-IOKR
%
% OUTPUTS:
% C:                regression model
%
%======================================================

    n_kx = length(KX_train_list); % number of input kernels
    n_train = size(KX_train_list{1},1);
    
    Ic = cell(n_kx,1);
    for k = 1:n_kx
        Ic{k} = eye(n_train);
    end
    I = cell2mat(Ic);
    clear Ic;
    
    KX_train = blkdiag(KX_train_list{1:n_kx});
    
    % Training MP-IOKR with the reverse IOKR
    A1 = I - M * Mean_KY_train_C;
    AAt = A1*A1' + M * Cov_KY_train_C * M;
    
    
    C = (lambda * eye(n_kx*n_train) + KX_train * AAt);

end



