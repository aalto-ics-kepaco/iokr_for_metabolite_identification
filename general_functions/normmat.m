function Kn = normmat( K, D1, D2 )
%======================================================
% DESCRIPTION:
% Normalization of a kernel Gram matrix
%
% INPUTS:
% K:    Gram matrix of a kernel k between two sets X1 and X2
% D1:   vector containing the diagonal of the Gram matrix of k between X1 and X1
% D2:   vector containing the diagonal of the Gram matrix of k between X2 and X2
%
% OUTPUT:
% Kn:   the normalized Gram matrix of k 
%
%======================================================

    if nargin == 1 
%         if isequal(K,K') % symmetric matrix (X1 = X2)
            D = diag(K);
            Kn = K ./ sqrt(D * D');
%         else
%             assert ('The function requires 3 arguments')
%         end
    else
        Kn = K ./ sqrt(D1 * D2');
    end

end

