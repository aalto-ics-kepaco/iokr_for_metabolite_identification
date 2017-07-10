function [ K_c ] = center( K, mean_K_train, center_opt, mean_K_1_train, mean_K_train_2 )
%======================================================
% DESCRIPTION:
% Centering a kernel Gram matrix in the feature space relative to the n
% feature vectors (training feature vectors in practice)
%
% INPUTS:
%  K:               Gram matrix of size n1*n2 between the sets X1 and X2
%  mean_K_train:    mean(K_train,2) where K_train is the Gram matrix over the training set
%  center:          binary value indicating if K has to be centered or not
%  mean_K_1_train:  Gram matrix of size n1*n_train between X1 and X_train
%  mean_K_train_2:  Gram matrix of size n_train*n2 between X_train and X2
%
% OUTPUT:
%  K_c:             centered Gram matrix of size n1 * n2
%
%======================================================

    [n1,n2] = size(K);
    
    if nargin == 3 % when X1 = X2 = X_train
        mean_K_1_train = mean_K_train';
        mean_K_train_2 = mean_K_train;
    end
    
    
    if center_opt == 1
        if all ([n1,n2] > 1) 
            % Kernel matrix provided
            K_c = K - repmat(mean_K_train_2,n1,1) - repmat(mean_K_1_train,1,n2) + mean(mean_K_train);
        else
            % Only diagonal of kernel matrix provided, e.g. if we need to
            % center the kernel matrix between all candidates.
            K_c = K(:) - mean_K_train_2(:) - mean_K_1_train(:) + mean(mean_K_train);
        end 
    else
        K_c = K;
    end
end


