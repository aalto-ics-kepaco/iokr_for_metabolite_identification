function [ gamma ] = select_gamma_entropy( X, KY_opt )
%======================================================
% DESCRIPTION:
% Find the gamma parameter of a Gaussian kernel such that the entropy of
% the Gaussian kernel is maximized
%
% INPUTS:
% X:                matrix of size d*n containing the input or output vectors
% base_kernel_type: type of the base kernel that will be used for computing the
%                   pairwise distances. Typically 'linear' or 'tanimoto'
%                   can be used, but not kernels that have parameter(s).
%                   The Gaussian kernel is then applied on top.
%
% OUTPUT:
% gamma:            Gaussian kernel parameter
%
%======================================================

    n = size(X,2);
    
    % computation of the base kernel (linear or Tanimoto)
    base_ker_opt = set_kernel_opt(KY_opt.base_kernel);
    K = build_kernel(X, X, base_ker_opt);
   
    % computation of all pairwise distances
    Kdist = diag(K)*ones(1,n) + ones(n,1)*diag(K)' - 2*K;

    % the gamma parameter maximising the entropy of the Gaussian kernel is
    % obtained by averaging all the pairwise distances and inverting this mean
    gamma = (n*(n-1)) / (2*sum(sum(triu(Kdist,1))));
        
end

