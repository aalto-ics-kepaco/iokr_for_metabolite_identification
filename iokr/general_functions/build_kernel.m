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
            
            K = build_linear_kernel( X1, X2 );
            
        case 'polynomial'
            
            K = build_polynomial_kernel(X1, X2, kernel_opt.scale, kernel_opt.offset, kernel_opt.degree);
            
        case 'tanimoto'
            
            K = build_tanimoto_kernel(X1, X2);

        case 'gaussian'
            
            switch kernel_opt.base_kernel
                case 'linear'
                    
                    Kb_1 = build_linear_kernel(X1, X1);
                    Kb_2 = build_linear_kernel(X2, X2);
                    Kb_12 = build_linear_kernel(X1, X2);

                case 'tanimoto'

                    Kb_1 = build_tanimoto_kernel(X1, X1);
                    Kb_2 = build_tanimoto_kernel(X2, X2);
                    Kb_12 = build_tanimoto_kernel(X1, X2);
            end
            
            % Gaussian kernel on top of the base kernel
            K = build_gaussian_kernel( Kb_1, Kb_2, Kb_12, kernel_opt.gamma);
            
    end

end

% linear kernel
function [ K_lin ] = build_linear_kernel( X1, X2 )
   
    K_lin = X1' * X2;

end

% polynomial kernel
function [ K_poly ] = build_polynomial_kernel( X1, X2, scale, offset, degree )
   
    K_poly = (scale * X1'*X2 + offset).^degree;

end

% Tanimoto kernel
function [ K_tan ] = build_tanimoto_kernel( X1, X2 )

    n1 = size(X1,2);
    n2 = size(X2,2);
    
    P = X1' * X2;
    D1 = sum(X1.*X1)';
    D2 = sum(X2.*X2)';

    K_tan = P ./ (D1 * ones(1,n2) + ones(n1,1) * D2' - P);

end

% Gaussian kernel
function [ K_gauss ] = build_gaussian_kernel( Kb_1, Kb_2, Kb_12, gamma)

    [n1,n2] = size(Kb_12);
    
    % Pairwise distances
    Dis = diag(Kb_1) * ones(1,n2) + ones(n1,1) * diag(Kb_2)' -2 * Kb_12;
            
    K_gauss = exp(-gamma * Dis);
    
end


