function [ vec_gamma_opt ] = Select_param_reverse_IOKR( KX_train_list, Psi_train, val_gamma )
%======================================================
% DESCRIPTION:
% Selection of the regularization parameters for the reverse IOKR approach 
% in the case of multiple input kernels using leave-one-out cross-validation
% (one parameter per input kernel is selected)
%
% INPUTS:
% KX_train_list:    list of input Gram matrices (size n_train*n_train)
% Psi_train:        matrix of size d*n_train containing the training output feature vectors
% val_gamma:        vector of stricly positive values among which the regularization 
%                   parameters of reverse IOKR will be selected
%
% OUTPUTS:
% vec_gamma_opt:    vector containing the different selected
%                   regularization parameters (one for each input kernel)
%
%======================================================

    n_kx = length(KX_train_list); % number of input kernels
    
    n_train = size(Psi_train,2);
    
    KY_train = Psi_train'*Psi_train;

    vec_gamma_opt = zeros(n_kx,1);
    for i = 1:n_kx
        
        mse = zeros(length(val_gamma),1);
        for ig = 1:length(val_gamma)
            gamma = val_gamma(ig);

            B = (gamma*eye(n_train) + KY_train) \ KY_train;
            LOOE = (eye(n_train)-B) / diag(diag(eye(n_train)-B));

            % Compute mean squared error
            mse(ig) = 1/n_train * trace(LOOE' * KX_train_list{i} * LOOE);
        end
        [~,ind_gamma_opt] = min(mse);
        vec_gamma_opt(i) = val_gamma(ind_gamma_opt);
        
    end

end

