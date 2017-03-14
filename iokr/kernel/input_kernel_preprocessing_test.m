function [ KX_train_test ] = input_kernel_preprocessing_test( KX_list_train_test, KX_list_test, train_process, ker_center )
%======================================================
% DESCRIPTION:
% Preprocessing of the input kernel matrices between the training set and
% the test set
%
% INPUTS:
% KX_list_train_test:   cell array of size: n_kx*1 containing the input kernel 
%                       matrices between the training and test sets
% KX_list_test:         cell array of size: n_kx*1 containing the input kernel 
%                       matrices on the test set
% train_process:        structure containing the information used for processing the 
%                       training samples
% ker_center:           binary value indicating if the kernel matrices should
%                       be centered or not
%
% OUTPUTS:
% KX_train_test:        linear combination of the preprocessed input kernel matrices 
%                       between the training and the test sets
%
%======================================================
    
    n_kx = length(KX_list_train_test); % number of input kernels
    
    [n_train, n_test] = size(KX_list_train_test{1}); % number of training and test examples
    
    KX_train_test = zeros(n_train, n_test);
    
    for i = 1:n_kx
        
        % Centering and normalization
        KX_train_test_i_cn = kernel_preprocessing_test(KX_list_train_test{i}, KX_list_test{i}, train_process(i), ker_center);
        
        % Kernel combination
        KX_train_test = KX_train_test + train_process(i).w * KX_train_test_i_cn;
    end


end

