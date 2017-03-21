function [ Psi_train, process_output ] = output_feature_preprocessing_train( Y_train, ker_center )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

    % Centering and normalization
    Psi_train = norma(Y_train, mean(Y_train,2), ker_center);
    
    process_output.mean = mean(Y_train,2);

end

