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
    

    n1 = size(X1,2);
    n2 = size(X2,2);
    
    P = X1' * X2; % inner product
            
    switch kernel_opt.type
        case 'linear'
            
            K = P;
            
        case 'tanimoto'
            
            D1 = sum(X1.*X1)';
            D2 = sum(X2.*X2)';
            
            K = P ./ (D1 * ones(1,n2) + ones(n1,1) * D2' - P);
            
        case 'tanimoto_gaussian'
            
            D1 = sum(X1.*X1)';
            D2 = sum(X2.*X2)';
            
            P1 = X1' * X1;
            P2 = X2' * X2;
            
            K_12 = P ./ (D1 * ones(1,n2) + ones(n1,1) * D2' - P);
            K_1 = P1 ./ (D1 * ones(1,n1) + ones(n1,1) * D1' - P1);
            K_2 = P2 ./ (D2 * ones(1,n2) + ones(n2,1) * D2' - P2);
            Dis = diag(K_1) * ones(1,n2) + ones(n1,1) * diag(K_2)' -2 * K_12;
            K = exp(-kernel_opt.gamma * Dis);
            
        case 'gaussian'
            
            D1 = sum(X1.*X1)';
            D2 = sum(X2.*X2)';
            
            Dis = D1 * ones(1,n2) + ones(n1,1) * D2' -2 * P;
            K = exp(- kernel_opt.gamma * Dis);
            
        case 'polynomial'
            
            K = (kernel_opt.scale * P + kernel_opt.offset).^kernel_opt.degree;

    end

end

