function [ M ] = Train_reverse_IOKR_kernel( KY_train, vec_gamma )
%======================================================
% DESCRIPTION:
% Training of reverse IOKR in the case of multiple input kernels and a
% kernel in output
%
% INPUTS:
% KY_train:     output Gram matrix on the training set (size n_train*n_train)
% gamma:        vector containing the different regularization parameters 
%               (one for each input kernel)
%
% OUTPUTS:
% M:            reverse IOKR model
%
%======================================================
    
    n_kx = length(vec_gamma); % number of input kernels
    
    n_train = size(KY_train,1);
    
    Mc = cell(n_kx,1);
    for k = 1:n_kx
        
        Mc{k} = inv(vec_gamma(k) * eye(n_train) + KY_train);
        
    end
    
    M = cell2mat(Mc);

end

