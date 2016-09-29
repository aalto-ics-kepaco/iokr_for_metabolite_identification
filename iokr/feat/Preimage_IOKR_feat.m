function [ scores ] = Preimage_IOKR_feat( Psi_pred, Y_C_test, mean_Y_train, ker_center )
%======================================================
% DESCRIPTION:
% Preimage of IOKR in the case of a feature represention in output
%
% INPUTS:
% Psi_pred:         matrix of size d*n_test containing the predicted output
%                   feature vectors of the test examples
% Y_C_test:         test candidate sets
% mean_Y_train:     mean of the training output vectors
% ker_center:       value of 0 of 1 indicating if the output feature/kernel
%                   should be centered.
%
% OUTPUTS:
% score:            cell of length n_test, the jth entry containing a vector of score 
%                   for the candidates of the test example j
%
%======================================================
    n_test = Y_C_test.getNumberOfExamples();
    
    % Pre-image
    scores = cell(n_test,1);
    for j = 1:n_test
        if (isnan (Y_C_test.getCandidateSet (j, 0, 'num')))
            scores{j} = NaN;
            
            continue;
        end % if
        
        Y_Cj = full (Y_C_test.getCandidateSet (j, 0, 'data'));
        Psi_Cj = norma(Y_Cj, mean_Y_train, ker_center); %  centering and normalization
        scores{j} = Psi_pred(:,j)' * Psi_Cj;         
    end
end