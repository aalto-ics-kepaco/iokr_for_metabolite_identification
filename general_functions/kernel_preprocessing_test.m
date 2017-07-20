function K_train_test_cn = kernel_preprocessing_test (K_train_test, K_test, ...
    train_process, ker_center, is_diagonal_K_test)
%======================================================
% DESCRIPTION:
% Preprocessing of a kernel matrix between training and test examples
%
% INPUTS:
% K_train_test:     kernel matrix between training and test examples
% K_test:           kernel matrix between test examples
% train_process:    structure containing the information needed for processing the 
%                   test samples identically to the training samples
% ker_center:       binary value indicating if the kernel matrices should
%                   be centered or not
% is_diagonal_K_test: Binary value indicating, whether the K_test variable
%                     contains only the diagonal elements of K_test, i.e. 
%                     K_test = diag(K_test).
% 
%
% OUTPUTS:
% K_train_test_cn:  preprocessed kernel matrix between the training and test exanples
%
%======================================================
    if (nargin < 5)
        is_diagonal_K_test = false;
    end % if

    % centering
    mean_K_train      = train_process.mean;
    mean_K_train_test = mean (K_train_test, 1);
    
    K_train_test_c = center (K_train_test, mean_K_train, ker_center, ...
        mean_K_train', mean_K_train_test);
    K_test_c       = center (K_test,       mean_K_train, ker_center, ...
        mean_K_train_test', mean_K_train_test, is_diagonal_K_test);
    
    % normalization    
    if (is_diagonal_K_test) 
        K_train_test_cn = normmat(K_train_test_c, train_process.diag_c, ...
            K_test_c);
    else
        K_train_test_cn = normmat(K_train_test_c, train_process.diag_c, ...
            diag (K_test_c));
    end % if
end

