function [ KX_train_test ] = mpiokr_input_kernel_preprocessing_test( KX_list_train_test, KX_list_test, process_input, mp_iokr_param, gamma_opt )
%======================================================
% DESCRIPTION:
% Preprocessing of the input kernel matrices between training and test
% examples in the case of MP-IOKR
%
% INPUTS:
% KX_list_train_test: cell array of size: n_kx*1 containing the input kernel 
%                     matrices between the training and test 
% KX_list_test:     cell array of size: n_kx*1 containing the test input kernel matrices
% process_input:    structure containing the information needed for processing the 
%                   test samples identically to the training samples
% mp_iokr_param:    1*1 struct array containing information relative to
%                   centering and multiple kernel learning
% gamma_opt:        selected value(s) for the regularization parameter gamma
%
% OUTPUTS:
% KX_train_test:    linear combination of the preprocessed input kernel matrices 
%                   between the training and test exanples
%
%======================================================
    
     switch mp_iokr_param.rev_iokr
            case 'joint' % joint approximation
                
                KX_train_test = {input_kernel_preprocessing_test(KX_list_train_test, KX_list_test, process_input, mp_iokr_param.center)};
                        
            case 'separate' % separate approximation
                
                switch nargin 
                    case 4
                        n_kx = length(KX_list_test);
                        KX_train_test = cell(n_kx, 1);
                        for j = 1:n_kx
                            KX_train_test{j} = input_kernel_preprocessing_test(KX_list_train_test{j}, KX_list_test{j}, process_input(j), mp_iokr_param.center);
                        end
                        
                    case 5
                          
                       [gamma_opt_u,~,ic] = unique(gamma_opt);

                        n_kx = length(gamma_opt_u);
                        KX_train_test = cell(n_kx, 1);
                        for j = 1:n_kx
                            KX_train_test{j} = input_kernel_preprocessing_test(KX_list_train_test(ic==j), KX_list_test(ic==j), process_input(j), mp_iokr_param.center);
                        end
                end
          
     end
end

