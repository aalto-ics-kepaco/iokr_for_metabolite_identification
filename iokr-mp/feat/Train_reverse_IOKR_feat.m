function [ M ] = Train_reverse_IOKR_feat( Psi_train, vec_gamma )
%======================================================
% DESCRIPTION:
% Training of reverse IOKR in the case of multiple input kernels
%
% INPUTS:
% Psi_train:    matrix of size d*n_train containing the training output feature vectors
% gamma:        vector containing the different regularization parameters 
%               (one for each input kernel)
%
% OUTPUTS:
% M:            reverse IOKR model
%
%======================================================

    n_kx = length(vec_gamma); % number of input kernels
    
    n_train = size(Psi_train,2);
    
    KY_train = Psi_train'*Psi_train;
    
    Mc = cell(n_kx,1);
    for k = 1:n_kx
        Mc{k} = (vec_gamma(k) * eye(n_train) + KY_train) \ (Psi_train');
    end
    
    M = cell2mat(Mc);
end

