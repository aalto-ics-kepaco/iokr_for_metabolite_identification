function [ KY_train_cand_cn ] = output_kernel_preprocessing_test( Y_train, Y_cand, KY_par, process_output, ker_center )
%======================================================
% DESCRIPTION:
% Preprocessing of the output kernel matrix between training and candidate
% examples
%
% INPUTS:
% Y_train:          matrix containing the training output vectors
% Y_train:          matrix containing the candidate output vectors
% KY_par:           structure containing the output kernel parameters
% process_output:   structure containing the information needed for processing the 
%                   test samples identically to the training samples
% ker_center:       binary value indicating if the kernel matrix should
%                   be centered or not
%
% OUTPUTS:
% KY_train_cand_cn: preprocessed output kernel matrix between training
%                   and candidate examples
% 
%
%======================================================

    KY_train_cand = build_kernel (Y_train, Y_cand, KY_par);
    KY_cand_diag  = build_kernel (Y_cand, Y_cand, KY_par, true);
    % Only the diagonal of the kernel matrix KY_cand is needed for
    % normalization.
    
    KY_train_cand_cn = kernel_preprocessing_test (KY_train_cand, KY_cand_diag, ...
        process_output, ker_center, ...
        true); % true ... only the diagonal matrix of KY_cand is provided.
    
    % Here we can put the varargout so that we return KY_cand_diag
end
