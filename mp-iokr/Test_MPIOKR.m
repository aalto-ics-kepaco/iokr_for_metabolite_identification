function [ score ] = Test_MPIOKR (KX_list_train_test, KX_list_test, train_model, ...
    Y_train, Y_C_train, Y_C_test, mp_iokr_param, ker_center, debug_param )
%======================================================
% DESCRIPTION:
% Prediction step of MP-IOKR
%
% INPUTS:
% KX_list_train_test:   cell array containing the input kernel matrices
%                       between the training and the test examples
% KX_list_test:         cell array containing the test input kernel matrices
% train_model:          1*1 struct array containing the regression model
%                       and information on the training data
% Y_train:              matrix of size d*n_train containing the training fingerprints
% Y_C_test:             cell array containing the candidate fingerprints
%                       (each element of the cell array corresponds to a candidate set)
% iokr_param:           1*1 struct array containing information relative to
%                       centering and multiple kernel learning
%
% OUTPUTS:
% score:                cell array containing the scores obtained for each
%                       candidate set
%
%======================================================
    KY_par_opt = train_model.KY_par;

    % Computation of the input kernel between training and test examples
    KX_list_train_test_processed = mpiokr_input_kernel_preprocessing_test ( ...
        KX_list_train_test, KX_list_test, train_model.process_input,        ...
        mp_iokr_param, train_model.gamma_opt);
    
    KX_train_test = cell2mat (KX_list_train_test_processed);
    
    % Prediction on the test set
    B = train_model.C * KX_train_test;
    
    % Pre-image
    if (~ strcmp (KY_par_opt.type, 'linear'))
        if (debug_param.verbose)
            tic;
        end % if;
        
        Y_train_all = [Y_train, ...
            cell2mat(arrayfun (@(idx) double (Y_C_train.getCandidateSet (idx, true, 'data')), ...
            1:Y_C_train.getNumberOfExamples(), 'UniformOutput', false))];
        n_train = size (Y_train, 2);
        n_Ci = arrayfun (@(idx) Y_C_train.getCandidateSet (idx, true, 'num'), ...
            1:Y_C_train.getNumberOfExamples());
        n_C_train = sum(n_Ci); % total number of candidates
        ind_S = 1:n_train;
        ind_C = n_train + (1:n_C_train);
        
        % Build the V and D matrices
        V = zeros(n_C_train, n_train);
        D = zeros(n_C_train, n_C_train);
        ind_0 = 0;
        for i = 1:n_train
            ind_i = ind_0+(1:n_Ci(i));
            V(ind_i,i) = 1;
            D(ind_i,ind_i) = 1/sqrt(n_Ci(i));
            ind_0 = ind_0+n_Ci(i);
        end
        V = sparse (V);
        D = sparse (D);
        
        D_squared              = D^2;
        D_I_minus_D_squaredVVt = D *(sparse (1:n_C_train, 1:n_C_train, 1)-D_squared*(V*V'))';
        
        KY_all_diag = build_kernel (Y_train_all, Y_train_all, KY_par_opt, true);        % O(l + m), only diagonal
        KY_S        = build_kernel (Y_train,     Y_train, KY_par_opt);                  % O(l^2)
        KY_all_S    = build_kernel (Y_train_all, Y_train, KY_par_opt);                  % O((l + m) * l)
        
        mean_KY_S     = mean (KY_S, 2);
        mean_KY_all_S = mean (KY_all_S, 2);
        clear KY_all_S KY_S;
        
        KY_all_diag_c = center (KY_all_diag, mean_KY_S, ker_center, ...
            mean_KY_all_S, mean_KY_all_S');
        clear KY_all_diag;
        
        if (debug_param.verbose)
            fprintf ('Scoring pre-processing: %.3fs\n', toc);
        end % if
    end % if
    
    n_test = Y_C_test.getNumberOfExamples();
    score  = cell (n_test,1);
    for j = 1:n_test
        if (debug_param.verbose)
            tic;
        end % if
        
        switch KY_par_opt.type
            case 'linear'
                                                 
                Psi_Cj = norma(Y_C_test.getCandidateSet (j, false, 'data'), mean(Y_train,2), ker_center);
                
                score{j} = B(:,j)' * Psi_Cj;
                
            otherwise                
                % Computation of the output kernel between the training +
                % candidate training examples and the candidate set Cj
                KY_Cj_diag  = build_kernel (Y_C_test.getCandidateSet (j, false, 'data'), ...    % O(m_j),   only diagonal
                    Y_C_test.getCandidateSet (j, false, 'data'), KY_par_opt, true);

                KY_all_Cj = build_kernel (Y_train_all, Y_C_test.getCandidateSet (j, false, 'data'), ...
                    KY_par_opt);                                                                % O((l + m) * m_j)
                
                % Centering  

                KY_Cj_diag_c  = center (KY_Cj_diag, mean_KY_S, ker_center, ...
                    mean (KY_all_Cj(ind_S, :), 1)' , mean (KY_all_Cj(ind_S, :), 1));
                
                KY_all_Cj_c = center (KY_all_Cj, mean_KY_S, ker_center, ...
                    mean_KY_all_S, mean (KY_all_Cj(ind_S, :), 1));
                
                clear KY_Cj_diag KY_all_Cj;
                
                % Normalization
                KY_all_Cj_cn = normmat (KY_all_Cj_c, KY_all_diag_c, KY_Cj_diag_c);
                clear KY_Cj_diag_c KY_all_Cj_c;
                
                % Tilde kernel 
                KY_tilde_Cj = [KY_all_Cj_cn(ind_S,:) - V'*D_squared*KY_all_Cj_cn(ind_C,:); ...
                    D_I_minus_D_squaredVVt*KY_all_Cj_cn(ind_C,:)];
                
                score{j} = B(:,j)' * KY_tilde_Cj;           
        end % switch
        
        if (debug_param.verbose)
           fprintf ('Scoring %d/%d with m=%d train-candidates and m_j=%d test-candidates: %.3fs\n', ...
               j, n_test, n_C_train, ...
               Y_C_test.getCandidateSet (j, false, 'num'), toc);
        end % if        
    end % for

end

