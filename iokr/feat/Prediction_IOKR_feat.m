function [ Psi_pred ] = Prediction_IOKR_feat( KX_train_test, C )
%======================================================
% DESCRIPTION:
% Test prediction of IOKR in the case of a feature represention in output
%
% INPUTS:
% KX_train_test:    input kernel matrix between training and test examples (size n_train*n_test)
% C:                IOKR model obtained with the function Train_IOKR_feat.m
%
% OUTPUTS:
% Psi_pred:         matrix of size d*n_test containing the predicted output
%                   feature vectors of the test examples
%
%======================================================
    Psi_pred = C * KX_train_test;
end

