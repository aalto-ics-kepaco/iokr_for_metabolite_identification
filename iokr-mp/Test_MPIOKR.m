function [ score ] = Test_MPIOKR( KX_list_train_test, KX_list_test, train_model, Y_train, Y_C_train, Y_C_test, mp_iokr_param )
%======================================================
% DESCRIPTION:
% Prediction step of IOKR
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

    ker_center = mp_iokr_param.center;

    KY_par_opt = train_model.KY_par;

    % Computation of the input kernel between training and test examples
    KX_list_train_test_processed = mpiokr_input_kernel_preprocessing_test(KX_list_train_test, KX_list_test, train_model.process_input, mp_iokr_param, train_model.gamma_opt);
    
    KX_train_test = cell2mat(KX_list_train_test_processed);
    
    % Prediction on the test set
    B  = train_model.C * KX_train_test;
    
    % Pre-image
    n_test = length(Y_C_test); % number of test examples
    score = cell(n_test,1);
    for j = 1:n_test
        
        switch KY_par_opt.type
            case 'linear'
                                                 
                Psi_Cj = norma(Y_C_test{j}, mean(Y_train,2), ker_center);
                
                score{j} = B(:,j)' * Psi_Cj;
                
            otherwise
                                
                Y_train_all = [Y_train,cell2mat(Y_C_train)];
                n_train = size(Y_train,2);
                n_Ci = cellfun(@(x) size(x,2),Y_C_train);
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
                
                % Computation of the output kernel between the training +
                % candidate training examples and the candidate set Cj
                KY_all = build_kernel(Y_train_all,Y_train_all, KY_par_opt);
                KY_all_Cj = build_kernel(Y_train_all,Y_C_test{j},KY_par_opt);
                KY_Cj = build_kernel(Y_C_test{j},Y_C_test{j},KY_par_opt);
                
                % Centering
                KY_all_c = center(KY_all, mean(KY_all(ind_S,ind_S),2), ker_center, mean(KY_all(:,ind_S),2), mean(KY_all(ind_S,:),1)); % centering
                KY_all_Cj_c = center(KY_all_Cj, mean(KY_all(ind_S,ind_S),2), ker_center, mean(KY_all(:,ind_S),2), mean(KY_all_Cj(ind_S,:),1) );
                KY_Cj_c = center(KY_Cj, mean(KY_all(ind_S,ind_S),2), ker_center, mean(KY_all_Cj(ind_S,:),1)' , mean(KY_all_Cj(ind_S,:),1));
                
                % Normalization
                KY_all_Cj_cn = normmat(KY_all_Cj_c,diag(KY_all_c),diag(KY_Cj_c));
                
                % Tilde kernel 
                KY_tilde_Cj = [KY_all_Cj_cn(ind_S,:) - V'*D^2*KY_all_Cj_cn(ind_C,:); D*(eye(n_C)-D^2*(V*V'))'*KY_all_Cj_cn(ind_C,:)];
                
                score{j} = B(:,j)' * KY_tilde_Cj;
        end
    end

end

