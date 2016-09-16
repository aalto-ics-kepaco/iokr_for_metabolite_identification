function [ score ] = Preimage_MP_IOKR_feat(Psi_pred, Y_C_test, mean_Y_train, ker_center)
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

    n_test = length(Y_C_test);
    
    % Pre-image
    score = cell(n_test,1);
    for j = 1:n_test
        Psi_Cj = norma(Y_C_test{j}, mean_Y_train, ker_center); %  centering and normalization
        score{j} = Psi_pred(:,j)' * Psi_Cj;         
    end

end