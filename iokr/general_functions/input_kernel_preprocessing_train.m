function [ KX_train, process_input ] = input_kernel_preprocessing_train( KX_list_train, w, ker_center )
%======================================================
% DESCRIPTION:
% Preprocessing of the training input kernel matrices
%
% INPUTS:
% KX_list_train:    cell array of size: n_kx*1 containing the input kernel 
%                   matrices on the training set
% w:                vector containing the mkl weights
% ker_center:       binary value indicating if the kernel matrices should
%                   be centered or not
%
% OUTPUTS:
% KX_train:         linear combination of the preprocessed input kernel matrices 
%                   on the training set
% process_input:    structure containing the information needed for processing the 
%                   test samples identically to the training samples
%
%======================================================
    
    n_kx = length(KX_list_train); % number of input kernels
    
    n_train = size(KX_list_train{1}); % number of training examples

    % Centering and normalization of the input kernels and kernel combination
    
    KX_train = zeros(n_train);
    process_input = struct('w',{},'mean',{},'diag_c',{});
    for i = 1:n_kx
        
        KX_i = KX_list_train{i};
        KX_i_c = center(KX_i, mean(KX_i,1), ker_center); % centering
        KX_i_cn = normmat(KX_i_c); % normalization
        
        KX_train = KX_train + w(i) * KX_i_cn;
        
        process_input(i).w = w(i);
        process_input(i).mean = mean(KX_i,1);
        process_input(i).diag_c = diag(KX_i_c);
        
    end
    
end

