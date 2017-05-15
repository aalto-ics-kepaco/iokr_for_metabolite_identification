function [ KY_train_cn, process_output ] = output_kernel_preprocessing_train( Y_train, KY_par, ker_center )
%======================================================
% DESCRIPTION:
% Preprocessing of the training output kernel matrix
%
% INPUTS:
% Y_train:          matrix containing the training output vectors
% KY_par:           structure containing the output kernel parameters
% ker_center:       binary value indicating if the kernel matrix should
%                   be centered or not
%
% OUTPUTS:
% KY_train_cn:      preprocessed output kernel matrix on the training set
% process_output:   structure containing the information needed for processing the 
%                   test samples identically to the training samples
%
%======================================================

    KY_train = build_kernel(Y_train, Y_train, KY_par);

    % Centering and normalization of the output kernel
    KY_train_c = center(KY_train, mean(KY_train,1), ker_center); % centering
    KY_train_cn = normmat(KY_train_c); % normalization
    
    process_output.mean = mean(KY_train,1);
    process_output.diag_c = diag(KY_train_c);
    
end
