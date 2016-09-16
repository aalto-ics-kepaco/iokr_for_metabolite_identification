function [ Mean_Psi_C_train, Cov_Psi_C_train ] = Compute_cov_mean_feat(Y_C_train, mean_Y, ker_center)
%======================================================
% DESCRIPTION:
% Computation of the candidate feature mean and covariance for
% magnitude-preserving IOKR
%
% INPUTS:
% Y_C_train:        training candidate sets
% mean_Y:           output vector mean of length d used for centering
%                   the candidate output vectors if ker_center is set to 1
% ker_center:       value of 0 of 1 indicating if the output feature/kernel
%                   should be centered.
%
% OUTPUTS:
% Mean_Psi_C_train: matrix of size d*n_train containing the averaged
%                   output feature vector of each training candidate set
% Cov_Psi_C_train:  matrix of size d*d containing the covariance of the
%                   candidate output feature vectors
%
%======================================================

    n_train = Y_C_train.getNumberOfExamples();
        
    d = size (Y_C_train.getCandidateSet (1, 0, 'data'), 1);
    Mean_Psi_C_train = zeros(d,n_train);
    Cov_Psi_C_train = zeros(d,d);
    for j = 1:n_train
        Y_Cj = Y_C_train.getCandidateSet(j, 1, 'data');
        nj = size(Y_Cj,2);

        % centering and normalization of the feature vectors
        Y_Cjn = norma(Y_Cj, mean_Y, ker_center);

        % computation of the mean and the covariance
        Mean_Psi_C_train(:,j) = mean(Y_Cjn,2);
        Cov_Psi_C_train = Cov_Psi_C_train + 1/nj*(Y_Cjn ...
            - repmat(Mean_Psi_C_train(:,j),1,nj))*(Y_Cjn - repmat(Mean_Psi_C_train(:,j),1,nj))';
    end

end


