function [ Psi_pred ] = Prediction_MP_IOKR_reverse_feat(KX_train_test_list, C)
%======================================================
% DESCRIPTION:
% Test prediction of MP-IOKR with reverse IOKR in the case of a feature
% represention in output
%
% INPUTS:
% KX_train_test:    list of input kernel matrices between training and test 
%                   examples (size n_train*n_test)
% C:                MP-IOKR model obtained with the function Train_MP_IOKR_reverse_feat.m
%
% OUTPUTS:
% Psi_pred:         matrix of size d*n_test containing the predicted output
%                   feature vectors of the test examples
%
%======================================================
    
    KX_train_test = cell2mat(KX_train_test_list);

    % Prediction on the test examples
    Psi_pred = C * KX_train_test;
end