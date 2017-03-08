function [ B ] = Prediction_IOKR_kernel( KX_train_test, C )
%======================================================
% DESCRIPTION:
% Test prediction of IOKR in the case of a kernel representation in output
%
% INPUTS:
% KX_train_test:    input kernel matrix between training and test examples (size n_train*n_test)
% C:      IOKR model obtained with the function Train_IOKR_kernel.m
%
% OUTPUTS:
% B:        Matrix used for the prediction of the output feature vectors
%
%======================================================

    B = C \ KX_train_test;

end

