function [ Mean_KY_train_C, Cov_KY_train_C ] = compute_cov_mean_kernel(Y_train, Y_C_train, OK)
%======================================================
% DESCRIPTION:
% Computation of the candidate kernel mean and covariance for
% magnitude-preserving IOKR
%
% INPUTS:
% Y_train:      matrix of size d*m containing the training output vectors
% Y_C_train:    cell of length n_train. Y_C_train{j} is a matrix of size d*n_j 
%               containing the output vectors of the candidates associated to y_j
% OK:           structure containing different informations about the
%               output kernel parameters and elements need for kernel centering
%   OK.center:  binary value indicating if the output kernel should be centered
%   OK.param:   output kernel parameter(s)
%   OK.mean:    mean of the training output kernel (before centering and normalization) used for centering
%   OK.diag:    diagonal of the training centered output kernel used for normalization
%
% OUTPUTS:
% Mean_KY_train_C: matrix of size n_train*n_train containing the averaged 
%                  output kernel values between the training set and the candidates
% Cov_KY_train_C:  matrix of size n_train*n_train containing the covariance of the 
%                  output kernel values between the training set and the candidates
%
%======================================================

    n_train = length(Y_C_train);
        
    Mean_KY_train_C = zeros(n_train,n_train);
    Cov_KY_train_C = zeros(n_train,n_train);
    for j = 1:n_train
        Y_Cj = Y_C_train{j};
        nj = size(Y_Cj,2);

        % Build the kernel between the training set and the candidate set
        KY_train_Cj = build_kernel_center_norm(Y_train, Y_Cj, OK.mean, OK.diag, OK.param, OK.center);

        % computation of the mean and the covariance
        Mean_KY_train_C(:,j) = mean(KY_train_Cj,2);
        KY_train_Cj_center = KY_train_Cj - repmat(mean(KY_train_Cj,2),1,nj);
        Cov_KY_train_C = Cov_KY_train_C + 1/nj* (KY_train_Cj_center*KY_train_Cj_center');
    end

end
