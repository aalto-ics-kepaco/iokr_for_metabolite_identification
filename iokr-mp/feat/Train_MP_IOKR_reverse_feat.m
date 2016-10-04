function [ C ] = Train_MP_IOKR_reverse_feat(KX_train_list, Psi_train, M, Mean_Psi_C_train, Cov_Psi_C_train, lambda)
%======================================================
% DESCRIPTION:
% Training of MP-IOKR with reverse IOKR in the case of a feature
% represention in output
%
% INPUTS:
% KX_train_list:        list containing one combined input kernel matrix on the training set (joint reverse IOKR) 
%                       or n_kx input kernel matrices on the training set, each of them being multiplied by the 
%                       corresponding mkl weight (separate reverse IOKR for each kernel)
% Psi_train:            matrix of size d*n_train containing the training output feature vectors
% M:                    list of 1 (joint case)/n_kx (separate case) reverse IOKR model(s)
% Mean_Psi_C_train:     matrix of size d*n_train containing the averaged 
%                       output feature vector of each training candidate set
% Cov_Psi_C_train:      matrix of size d*d containing the covariance of the 
%                       candidate output feature vectors
% lambda:               regularization parameter of MP-IOKR (>0)
%
% OUTPUTS:
% C:                    regression model
%
%======================================================

    n_kx = length(KX_train_list); % number of input kernels
    n_train = size(Psi_train,2);
    
    Ic = cell(n_kx,1);
    for k = 1:n_kx
        Ic{k} = eye(n_train);
    end
    I = cell2mat(Ic);
    clear Ic;
    
    KX_train = blkdiag(KX_train_list{1:n_kx});

    % Training MP-IOKR with the reverse IOKR
    A1 = I - M * Mean_Psi_C_train;
    PsiAt = (Psi_train - Mean_Psi_C_train) * A1' + Cov_Psi_C_train*M';
    AAt = A1 * A1' + M * Cov_Psi_C_train * M';
    
    C = PsiAt / (lambda * eye(n_kx*n_train) + KX_train * AAt);

end