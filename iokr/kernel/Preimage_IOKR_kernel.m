function [ score ] = Preimage_IOKR_kernel( Y_train, KY_par, B, Y_C_test, param )
%======================================================
% DESCRIPTION:
% Preimage of IOKR in the case of a kernel represention in output
%
% INPUTS:
% Y_train:          matrix of size d*n_train containing the training output
%                   vectors
% KY_par:           parameters of the output kernel
% B:                matrix equal to (lambda*eye(n_train) + KX_train) \ KX_train_test
% Y_C_test:         test candidate sets
% param:            structure containing the MP-IOKR parameters
%   param.center:       binary value indicating if the input and output
%                       kernel/feature vectors should be centered (1) or not (0)
%   param.mkl:          string indicating the MKL algorithm used for kernel combination 
%                       ('alignf' or 'unimkl')
%   param.cv:           string indicating the type of cross-validation
%                       ('cv' or 'loocv') for parameter selection
%
% OUTPUTS:
% score:            cell of length n_test, the jth entry containing a vector of score 
%                   for the candidates of the test example j
%
%======================================================

    n_test = length(Y_C_test);

    % Computation of the output kernel
    KY_train = build_kernel(Y_train, Y_train, KY_par);
    mean_KY_train = mean(KY_train,1);
    KY_train_c = center(KY_train, mean_KY_train, param.center); % centering and normalization
    
    
    % Pre-image
    score = cell(n_test,1);
    for j = 1:n_test        
        KY_train_Cj = build_kernel_center_norm(Y_train, Y_C_test{j}, ...
                    mean_KY_train, diag(KY_train_c), KY_par, param.center);
        score{j} = B(:,j)' * KY_train_Cj;
    end


end

