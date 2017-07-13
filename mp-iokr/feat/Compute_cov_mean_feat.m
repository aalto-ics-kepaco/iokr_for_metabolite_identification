function [ Mean_Psi_C_train, Cov_Psi_C_train ] = Compute_cov_mean_feat ( ...
    Y_C_train, mean_Y, ker_center, verbose)
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
% verbose:          binary indicating whether the progress of the
%                   calculation of the statistics should be printed
%
% OUTPUTS:
% Mean_Psi_C_train: matrix of size d*n_train containing the averaged
%                   output feature vector of each training candidate set
% Cov_Psi_C_train:  matrix of size d*d containing the covariance of the
%                   candidate output feature vectors
%
%======================================================
    if (nargin < 4)
        verbose = false;
    end % if
    
    if (verbose)
        reverseStr = '';
    end % if

    n_train = Y_C_train.getNumberOfExamples();
        
    d = size (Y_C_train.getCandidateSet (1, 0, 'data'), 1);
    Mean_Psi_C_train = zeros(d,n_train);
    Cov_Psi_C_train = zeros(d,d);
    for idx = 1:n_train
         if (verbose)
            percentDone = 100 * idx / n_train;
            msg = sprintf ('Calculate candidate statistics statistics: %d/%d (%3.1f perc)', idx, n_train, percentDone);
            fprintf ([reverseStr, msg]);
            reverseStr = repmat (sprintf ('\b'), 1, length (msg));           
        end % if
        
        nj = Y_C_train.getCandidateSet (idx, true, 'num');
        if (isnan (nj))
            % No candidate set for the desired example available.            
            continue
        end
        if (nj == 0)
            % If a candidate set does not contain any candidate no mean and
            % no covariance matrice needs to be calculated.
            continue;
        end % if
        
        Y_Cj = full (Y_C_train.getCandidateSet (idx, true, 'data'));       
        d = size (Y_Cj, 1);
        
        assert (d == numel (mean_Y), ...
            'Dimension of the "data" does not match the dimension of the mean vector.')
        assert (size (Y_Cj, 2) == nj);
        
        % centering and normalization of the feature vectors
        Psi_Cjn = norma(Y_Cj, mean_Y, ker_center);
        clear Y_Cj;

        % computation of the mean and the covariance
        Mean_Psi_C_train(:,idx) = mean(Psi_Cjn,2);
        Cov_Psi_C_train = Cov_Psi_C_train + 1/nj*(Psi_Cjn ...
            - repmat(Mean_Psi_C_train(:,idx),1,nj))*(Psi_Cjn - repmat(Mean_Psi_C_train(:,idx),1,nj))';
    end
    
    if (verbose)
        fprintf ('\n');
    end % if
end


