function [ w ] = mkl_weight( mkl_type, KX_list, KY )
%======================================================
% DESCRIPTION:
% Multiple Kernel Learning implementation for Alignf and Unimkl
%
% INPUTS:
% mkl_type: string indicating the MKL algorithm used for kernel combination 
%           ('alignf' or 'unimkl')
% KX_list:  cell array containing n_kx input Gram matrices (each of size n*n)
% KY:       Gram matrix of the output kernel (size: n*n)
%
% OUTPUTS:
% w:        vector containing the learnt weights (length: n_kx)
%
%======================================================

    n_kx = length(KX_list); % number of input kernels
    
    % Weight learning
    switch mkl_type
        case 'alignf'
            
            % Centering input and output training Gram matrices                           
            KX_c = cellfun(@(x) center(x,mean(x,1),1) ,KX_list, 'UniformOutput', false); 
            KY_c = center(KY, mean(KY,1),1);

            % Computation of centered alignments
            M = zeros(n_kx,n_kx);
            for i = 1:n_kx
                for j = i:n_kx
                    M(i,j) = sum(sum(KX_c{i}.*KX_c{j}));
                    M(j,i) = M(i,j);
                end
            end

            a = cellfun(@(x) sum(sum(x.*KY_c)), KX_c, 'UniformOutput', true);

            % Learning the weights w
            opts = optimoptions('quadprog','Display','off');
            w = quadprog(M,-a,[],[],[],[],zeros(n_kx,1),[],[],opts);
            w = w ./ norm(w);
    
        case 'unimkl'
            w = ones(n_kx,1)./n_kx;
    end

end

