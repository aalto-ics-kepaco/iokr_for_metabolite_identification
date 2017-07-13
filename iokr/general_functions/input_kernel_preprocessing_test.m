function [ KX_train_test ] = input_kernel_preprocessing_test( KX_list_train_test, KX_list_test, process_input, ker_center )
%======================================================
% DESCRIPTION:
% Preprocessing of the input kernel matrices between training and test examples
%
% INPUTS:
% KX_list_train_test: cell array of size: n_kx*1 containing the input kernel 
%                     matrices between the training and test 
% KX_list_test:     cell array of size: n_kx*1 containing the test input kernel matrices
% process_input:    structure containing the information needed for processing the 
%                   test samples identically to the training samples
% ker_center:       binary value indicating if the kernel matrices should
%                   be centered or not
%
% OUTPUTS:
% KX_train_test:    linear combination of the preprocessed input kernel matrices 
%                   between the training and test exanples
%
%======================================================

    KX_train_test = zeros(size(KX_list_train_test{1}));
    
    for i = 1:length(KX_list_train_test)
        % Centering and normalization
        KX_train_test_i_cn = kernel_preprocessing_test(KX_list_train_test{i}, KX_list_test{i}, process_input(i), ker_center);
        
        KX_train_test = KX_train_test + process_input(i).w * KX_train_test_i_cn;
    end
    
end

