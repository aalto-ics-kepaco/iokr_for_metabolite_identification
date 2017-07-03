function [ KX_train, process_input ] = mpiokr_input_kernel_preprocessing_train( KX_list_train, w, mp_iokr_param, gamma_opt )
%======================================================
% DESCRIPTION:
% Preprocessing of the input kernel matrices between training
% examples in the case of MP-IOKR
%
% INPUTS:
% KX_list_train:    cell array of size: n_kx*1 containing the training input kernel matrices
% w:                MKL weights
% mp_iokr_param:    1*1 struct array containing information relative to
%                   centering and multiple kernel learning
% gamma_opt:        selected value(s) for the regularization parameter gamma
%
% OUTPUTS:
% KX_train:         linear combination of the preprocessed input kernel matrices 
%                   between the training exanples
% process_input:    structure containing the information needed for processing the 
%                   test samples identically to the training samples
%
%======================================================

    switch mp_iokr_param.rev_iokr
        case 'joint' % joint approximation

            [KX_train_combined, process_input] = input_kernel_preprocessing_train(KX_list_train, w, mp_iokr_param.center);
            KX_train = {KX_train_combined};

        case 'separate' % separate approximation
            
            switch nargin
                case 3
                    
                    n_kx = length(KX_list_train);
                    KX_train = cell(n_kx, 1);
                    process_input = struct('w',{},'mean',{},'diag_c',{});
                    for j = 1:n_kx
                        [KX_train{j},process_input(j)] = input_kernel_preprocessing_train(KX_list_train{j}, w(j), mp_iokr_param.center);
                    end
                    
                case 4 % if gamma_opt is an argument of the function, then we group the kernels by unique gamma values
            
                   [gamma_opt_u,~,ic] = unique(gamma_opt);

                    n_kx = length(gamma_opt_u);
                    KX_train = cell(n_kx, 1);
                    process_input = struct('w',{},'mean',{},'diag_c',{});
                    for j = 1:n_kx
                        [KX_train{j},process_input(j)] = input_kernel_preprocessing_train(KX_list_train(ic==j), w(ic==j), mp_iokr_param.center);
                    end
            end
            
    end

end

