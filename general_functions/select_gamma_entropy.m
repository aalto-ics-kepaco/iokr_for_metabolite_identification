function [ gamma ] = select_gamma_entropy( X, kernel_opt )
%======================================================
% DESCRIPTION:
% Find the gamma parameter of a Gaussian kernel such that the entropy of
% the Gaussian kernel is maximized
%
% INPUTS:
% X:            matrix of size d*n containing the input or output vectors
% kernel_opt:   structure containing the kernel type and parameter(s)
%               (the kernel type should be Gaussian)
%
% OUTPUT:
% gamma:        Gaussian kernel parameter
%
%======================================================

    n = size(X,2);
    
    if strcmp(kernel_opt.type, 'gaussian')
        % computation of the base kernel (linear or Tanimoto)
        base_ker_opt = set_kernel_opt(kernel_opt.base_kernel);

        K = build_kernel(X, X, base_ker_opt);

        % computation of all pairwise distances
        Kdist = diag(K)*ones(1,n) + ones(n,1)*diag(K)' - 2*K;

        % the gamma parameter maximising the entropy of the Gaussian kernel is
        % obtained by averaging all the pairwise distances and inverting this mean
        gamma = (n*(n-1)) / (2*sum(sum(triu(Kdist,1))));
    else
        error ('select_gamma_entropy:InvalidArguments', ...
            'Entropy measure is only useful if the kernel type is gaussian: %s', ...
            kernel_opt.type);
    end % if
end % function

