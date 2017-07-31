function [ K ] = build_kernel( X1, X2, kernel_opt, only_diag )
%======================================================
% DESCRIPTION:
% Compute the gram matrix between the sets X1 and X2
%
% INPUTS:
% X1:           matrix of size d * n1
% X2:           matrix of size d * n2
% kernel_opt:   structure containing the kernel type and parameter(s). It can
%               be defined using the set_kernel_opt.m function
% only_diag:    only the diagonal elements of the kernel matrix are
%               calculated. This option can only be used when X1 = X2.
% 
% OUTPUTS:
% K:           Gram matrix of size n1 * n2
%
%======================================================

    if nargin < 4
        only_diag = false;
    end % if
    
    if only_diag && (~ isequal (X1, X2))
        error ('Diagonal of the kernel matrix makes only sense if both matrices are identical.');
    end % if

    % FIXME: Problems with logical matrices as '*' not defined.
    X1 = double (X1);
    X2 = double (X2);

    switch kernel_opt.type
        case 'linear'
            
            K = build_linear_kernel( X1, X2, only_diag );
            
        case 'polynomial'
            
            K = build_polynomial_kernel(X1, X2, ...
                kernel_opt.scale, kernel_opt.offset, kernel_opt.degree, ...
                only_diag);
            
        case 'tanimoto'
            
            K = build_tanimoto_kernel( X1, X2, only_diag );

        case 'gaussian'
            
            if only_diag
                K = ones (size (X1, 2), 1);
                return;
            end % if 
            
            switch kernel_opt.base_kernel
                case 'linear'

                    Kb_1_diag = build_linear_kernel(X1, X1, true);
                    Kb_2_diag = build_linear_kernel(X2, X2, true);
                    Kb_12     = build_linear_kernel(X1, X2, false);

                case 'tanimoto'

                    Kb_1_diag = build_tanimoto_kernel(X1, X1, true);
                    Kb_2_diag = build_tanimoto_kernel(X2, X2, true);
                    Kb_12     = build_tanimoto_kernel(X1, X2, false);
            end

            % Gaussian kernel on top of the base kernel
            K = build_gaussian_kernel( Kb_1_diag, Kb_2_diag, Kb_12, ...
                kernel_opt.gamma);

    end

end

% linear kernel
function [ K_lin ] = build_linear_kernel( X1, X2, only_diag)

    if only_diag
        K_lin = sum (X1 .* X2)';
    else
        K_lin = X1' * X2;
    end

end

% polynomial kernel
function [ K_poly ] = build_polynomial_kernel( X1, X2, scale, offset, degree, only_diag )
   
    if only_diag
        K_poly = (scale * sum (X1 .* X2)' + offset).^degree;
    else
        K_poly = (scale * X1'*X2 + offset).^degree;
    end

end

% Tanimoto kernel
function [ K_tan ] = build_tanimoto_kernel( X1, X2, only_diag )

    if only_diag
        K_tan = ones (size (X1, 2), 1);
        return;
    end % if 
    
    n1 = size(X1,2);
    n2 = size(X2,2);

    P = X1' * X2;
    D1 = sum(X1.*X1)';
    D2 = sum(X2.*X2)';

    K_tan = P ./ (D1 * ones(1,n2) + ones(n1,1) * D2' - P);

end

% Gaussian kernel
function [ K_gauss ] = build_gaussian_kernel( Kb_1_diag, Kb_2_diag, ...
    Kb_12, gamma )

    [n1,n2] = size(Kb_12);
    
    % Pairwise distances
    Dis = Kb_1_diag * ones(1,n2) + ones(n1,1) * Kb_2_diag' -2 * Kb_12;
            
    K_gauss = exp(-gamma * Dis);
    
end


