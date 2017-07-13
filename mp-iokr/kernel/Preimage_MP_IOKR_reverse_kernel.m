function [ score ] = Preimage_MP_IOKR_reverse_kernel(M, Y_train, OK, Y_C_train, Y_C_test, B, Mean_KY_train_C)
%======================================================
% DESCRIPTION:
% Test prediction of MP-IOKR with reverse IOKR in the case of a kernel
% represention in output
%
% INPUTS:
% M:
% Y_train:          matrix of size d*n_train containing the training outputfeature vectors
% OK:               structure containing different informations about the
%                   output kernel parameters and elements need for kernel centering
%   OK.center:      binary value indicating if the output kernel should be centered
%   OK.param:       output kernel parameter(s)
%   OK.mean:        mean of the training output kernel (before centering and normalization) used for centering
%   OK.diag:        diagonal of the training centered output kernel used for normalization
% Y_C_train:        training candidate sets
% Y_C_test:         test candidate sets 
% B:                MP-IOKR model obtained with the function Train_MP_IOKR_reverse_kernel.m
% Mean_KY_train_C:  matrix of size n_train*n_train containing the averaged 
%                   output kernel values between the training set and the candidates
%
% OUTPUTS:
% score:            cell of length n_test, the jth entry containing a vector of score 
%                   for the candidates of the test example j
%
%======================================================

    n_train = length(Y_C_train);
    n_test = length(Y_C_test);
    
    Ic = cell(n_kx,1);
    for k = 1:n_kx
        Ic{k} = eye(n_train);
    end
    I = cell2mat(Ic);
    clear Ic;
        
    % Pre-image
    score = cell(n_test,1);
    for j = 1:n_test 
        nj = size(Y_C_test{j},2);
        
        KY_train_Cj = build_kernel_center_norm(Y_train, Y_C_test{j}, OK.mean, OK.diag, OK.param, OK.center);
        mean_KY_C_Cj = zeros(n_train, nj);
        cov = zeros(n_train, nj);
        for i = 1:n_train
            ni = size(Y_train{i},2);
            KY_train_Ci = build_kernel_center_norm(Y_train, Y_C_train{i}, OK.mean, OK.diag, OK.param, OK.center);
            KY_Ci_Cj = build_kernel_center_norm(Y_C_train{i}, Y_C_test{j}, OK.mean, OK.diag, OK.param, OK.center);
            mean_KY_C_Cj(i,:) = mean(KY_Ci_Cj,1);
            cov = cov + 1/ni * (KY_train_Ci - repmat(Mean_KY_train_C(:,i),1,ni)) * (KY_Ci_Cj - repmat(mean_KY_C_Cj(i,:),1,ni));
        end
        F = (I - M*Mean_KY_train_C)*(KY_train_Cj - mean_KY_C_Cj) + M * cov;
        
        score{j} = B(:,j)' * F;
    end

end

