function [ K ] = build_kernel( X1, X2, kernel_opt )
%======================================================
% DESCRIPTION:
% Compute the gram matrix between the sets X1 and X2
%
% INPUTS:
% X1:           matrix of size d * n1
% X2:           matrix of size d * n2
% kernel_opt:   structure containing the kernel type and parameter(s). It can
%               be defined using the set_kernel_opt.m function
% 
% OUTPUTS:
% K:           Gram matrix of size n1 * n2
%
%======================================================
    
    switch kernel_opt.type
        case 'linear'
            
            K = X1' * X2;
            
        case 'gaussian'
            n1 = size(X1,2);
            n2 = size(X2,2);
            
            D1 = sum(X1.*X1)';
            D2 = sum(X2.*X2)';
            
            Dis = D1 * ones(1,n2) + ones(n1,1) * D2' -2 * X1' * X2;
            K = exp(- kernel_opt.gamma * Dis);
            
        case 'polynomial'
            
            K = (kernel_opt.scale * (X1' * X2) + kernel_opt.offset).^kernel_opt.degree;

    end

end

