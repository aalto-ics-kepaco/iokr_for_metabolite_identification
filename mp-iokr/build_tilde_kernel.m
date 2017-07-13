function [ KY_tilde, KY_S_C, V, D ] = build_tilde_kernel( Y, Y_C_list, KY_par, mp_iokr_param )
%======================================================
% DESCRIPTION:
% MP-IOKR: build the tilde output kernel matrix
%
% INPUTS:
% Y:            matrix of size d*n containing the training output vectors
% Y_C_list:     training candidate output vectors
% KY_par:       output kernel parameters
% mp_iokr_param:1*1 struct array containing information relative to
%               centering and multiple kernel learning
%
% OUTPUT:
% KY_tilde:     tilde output kernel matrix
% KY_S_C:       output kernel matrix between the training exmples and the
%               training candidate examples
% V:            matrix V of size n_C * n which indicates the correspondance between
%               the training examples and the training candidate examples
% D:            block diagonal matrix of size n_C*n_C
%
%======================================================
    
    n = size(Y,2);
    
    n_Ci = cellfun(@(x) size(x,2),Y_C_list);
    n_C = sum(n_Ci); % total number of candidates
    
    % Build the V and D matrices
    V = zeros(n_C, n);
    D = zeros(n_C, n_C);
    ind_0 = 0;
    for i = 1:n
        ind_i = ind_0+(1:n_Ci(i));
        V(ind_i,i) = 1;
        D(ind_i,ind_i) = 1/sqrt(n_Ci(i));
        ind_0 = ind_0+n_Ci(i);
    end
    
    ind_S = 1:n;
    ind_C = n + (1:n_C);
    
    Y_C = cell2mat(Y_C_list);
    
    % Building and processing the output kernel matrix
    KY1 = build_kernel([Y, Y_C], [Y, Y_C], KY_par);
    KY_c = center(KY1, mean(KY1(ind_S,ind_S),2), mp_iokr_param.center, mean(KY1(:,ind_S),2), mean(KY1(ind_S,:),1)); % centering
    KY = normmat(KY_c); % normalization
    
    KY_S_C = KY(ind_S,ind_C);
    
    KY_tilde = zeros(size(KY));
    KY_tilde(ind_S,ind_S) = KY(ind_S,ind_S) - V'*D^2*KY(ind_C,ind_S) - KY(ind_S,ind_C)*D^2*V + V'*D^2*KY(ind_C,ind_C)*D^2*V;
    KY_tilde(ind_S,ind_C) = (KY(ind_S,ind_C) - V'*D^2*KY(ind_C,ind_C)) * (eye(n_C) - D^2*(V*V'))*D;
    KY_tilde(ind_C,ind_S) = KY_tilde(ind_S,ind_C)';
    KY_tilde(ind_C,ind_C) = D*(eye(n_C) - D^2*(V*V'))'*KY(ind_C,ind_C)*(eye(n_C) - D^2*(V*V'))*D;

end

