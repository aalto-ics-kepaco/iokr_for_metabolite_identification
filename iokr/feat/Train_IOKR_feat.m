function [ C ] = Train_IOKR_feat( KX_train, Psi_train, lambda )
%======================================================
% DESCRIPTION:
% Training of IOKR in the case of a feature represention in output
%
% INPUTS:
% KX_train:     input kernel matrix on the training set (size n_train*n_train)
% Psi_train:    matrix of size d*n_train containing the training output feature vectors
% lambda:       regularization parameter (>0)
%
% OUTPUTS:
% C:            regression model
%
%======================================================

    n_train = size(KX_train,1);
    
    C = Psi_train / (lambda*eye(n_train) + KX_train);

end

