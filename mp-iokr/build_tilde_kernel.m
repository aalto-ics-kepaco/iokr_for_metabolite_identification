function [ KY_tilde, KY_S_C, V, D ] = build_tilde_kernel( Y, Y_C, KY_par, ...
    mp_iokr_param, debug_param)
    %======================================================
    % DESCRIPTION:
    % MP-IOKR: build the tilde output kernel matrix
    %
    % INPUTS:
    % Y:            matrix of size d*n containing the training output vectors
    % Y_C:          training candidates encapsulated in a CandidateSet object
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
    
    if (debug_param.verbose)
        tic;
    end % if

    n = size(Y,2);
    
    n_Ci = arrayfun (@(id) Y_C.getCandidateSet (id, true, 'num'), ...
        1:Y_C.getNumberOfExamples());
    % If a training example does not have a candidate we set its number of
    % candidates to zero. 
    n_Ci(isnan (n_Ci)) = 0; 
    n_C = sum(n_Ci); 
    
    % Build the V and D matrices
    V = zeros (n_C, n);
    ind_0 = 0;
    for i = 1:n
        ind_i = ind_0+(1:n_Ci(i));
        V(ind_i,i) = 1;
        ind_0 = ind_0+n_Ci(i);
    end
    V = sparse (V);
    D = sparse (diag (rude (n_Ci, 1 ./ sqrt (n_Ci))));
    
    ind_S = 1:n;
    ind_C = n + (1:n_C);
    
    
    % FIXME: Seems that not all candidate sets are logical.
    Y_C_data = arrayfun (@(id) double (Y_C.getCandidateSet (id, true, 'data')), ...
        1:Y_C.getNumberOfExamples(), 'UniformOutput', false);
    % If a training example does not have a candidate we do not consider
    % its feature vectors.
    has_Y_C  = cellfun (@(x) any (~ isnan (x(:))), Y_C_data);
    Y_C_data = cell2mat (Y_C_data(has_Y_C));
    
    % Building and processing the output kernel matrix
    KY1  = build_kernel ([Y, Y_C_data], [Y, Y_C_data], KY_par); % O((l + m)^2)
    KY_c = center (full (KY1), mean(KY1(ind_S,ind_S),2), mp_iokr_param.center, ...
        mean(KY1(:,ind_S),2), mean(KY1(ind_S,:),1)); % centering
    clear KY1; 
    KY   = normmat (KY_c); % normalization
    clear KY_c;
    
    KY_tilde = zeros (size(KY));
    % Some precalculations
    D_squared              = D^2;
    VVt                    = V * V';
    VtD_squared            = V' * D_squared;
    VtD_squaredK_CS        = VtD_squared * KY(ind_C,ind_S);
    VtD_squaredK_CC        = VtD_squared * KY(ind_C,ind_C);
    D_squaredVVt           = D_squared * VVt;
    I_minus_D_squaredVVt_D = (sparse (1:n_C, 1:n_C, 1) - D_squaredVVt) * D;

    KY_tilde(ind_S,ind_S) = KY(ind_S,ind_S) - VtD_squaredK_CS - ...
        VtD_squaredK_CS' + VtD_squaredK_CC*VtD_squared';
    clear VtD_squaredK_CS VtD_squared;

    KY_tilde(ind_S,ind_C) = (KY(ind_S,ind_C) - VtD_squaredK_CC) * I_minus_D_squaredVVt_D;
    clear VtD_squaredK_CC;

    KY_tilde(ind_C,ind_S) = KY_tilde(ind_S,ind_C)';
    KY_tilde(ind_C,ind_C) = I_minus_D_squaredVVt_D'*KY(ind_C,ind_C)*I_minus_D_squaredVVt_D;
    clear VVt D_squared D_squaredVVt I_minus_D_squaredVVt_D;
    
    if (debug_param.verbose)
        fprintf ('build_tilde_kernel (n = %d, n_C = %d): %.3fs\n', n, n_C, toc);
    end % if
    
    KY_S_C = KY(ind_S,ind_C);
end

