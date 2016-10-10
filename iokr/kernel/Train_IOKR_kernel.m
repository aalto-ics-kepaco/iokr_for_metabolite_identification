function [ C ] = Train_IOKR_kernel( KX_train, lambda )
%======================================================
% DESCRIPTION:
% Training of IOKR in the case of a kernel represention in output
%
% INPUTS:
% KX_train:     input kernel matrix on the training set (size n_train*n_train)
% lambda:       regularization parameter (>0)
%
% OUTPUTS:
% C:            (partial) regression model
%
%======================================================

    n_train = size(KX_train,1);
    
    C = lambda*eye(n_train) + KX_train;
    
end

