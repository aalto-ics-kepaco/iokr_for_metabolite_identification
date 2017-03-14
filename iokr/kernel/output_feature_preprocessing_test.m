function [ Psi_cand ] = output_feature_preprocessing_test( Y_cand, train_process, ker_center )
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

    Psi_cand = norma(Y_cand, train_process.mean, ker_center);

end

