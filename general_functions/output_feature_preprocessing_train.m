function [ Psi_train, process_output ] = output_feature_preprocessing_train( Y_train, ker_center )
%======================================================
% DESCRIPTION:
% Preprocessing of the training output feature vectors
%
% INPUTS:
% Y_train:          matrix containing the training output vectors
% ker_center:       binary value indicating if the feature vectors should
%                   be centered or not
%
% OUTPUTS:
% Psi_train:        preprocessed training output feature vectors
% process_output:   structure containing the information needed for processing the 
%                   test samples identically to the training samples
%
%======================================================

    % Centering and normalization
    Psi_train = norma(Y_train, mean(Y_train,2), ker_center);
    
    process_output.mean = mean(Y_train,2);

end

