function [ scores ] = Preimage_MP_IOKR_feat(Psi_pred, Y_C_test, mean_Y_train, ker_center)
%======================================================
% DESCRIPTION:
% Preimage of MP-IOKR with reverse IOKR in the case of a feature
% represention in output
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
        Y_Cj = Y_C_test.getCandidateSet (j, 0, 'data');
        if (any (any (isnan (Y_Cj))))
            % No candidate set for the desired example available. 
            warning ('No candidate set for desired example.');
            
            scores{j} = NaN;
            
            continue;
        end % if
             
        Psi_Cj = norma(Y_Cj, mean_Y_train, ker_center); %  centering and normalization
        scores{j} = Psi_pred(:,j)' * Psi_Cj;         
    end
end