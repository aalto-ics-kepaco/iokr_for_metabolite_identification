function [ scores, process_output, varargout ] = Test_IOKR( KX_list_train_test, KX_list_test, ...
    train_model, Y_train, Y_C_test, ker_center )
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
% Y_train:              matrix of size d*n_train containing the training output vectors
% Y_C_test:             cell array containing the candidate output vectors
%                       (each element of the cell array corresponds to a candidate set)
% ker_center:           binary value indicating if the input/output kernels should be
%                       centered or not
%
% OUTPUTS:
% score:                cell array containing the scores obtained for each
%                       candidate set
%
%======================================================
    calculate_squared_norm_of_prediction = (nargout == 3);

    % Computation of the input kernel between training and test examples
    KX_train_test = input_kernel_preprocessing_test (KX_list_train_test, ...
        KX_list_test, train_model.process_input, ker_center);
    
    % Prediction on the test set 
    % B = (lambda_opt * I + KX_train)^(-1) * KX_train_test
    switch train_model.model_representation
        case 'only_C'
            B = train_model.C \ KX_train_test;
        case 'Chol_decomp_of_C'
            y = linsolve (train_model.C,  KX_train_test, struct ('LT', true));
            B = linsolve (train_model.C', y,             struct ('UT', true));
    end % switch
    
    % Pre-image
    
    % Preprocessing of the training outputs    
    switch train_model.KY_par.representation
        case 'feature'
            switch train_model.KY_par.type 
                case 'linear'
                    % Y_train = Y_train
                case 'gaussian'
                    Y_train = train_model.KY_par.rff.getRandomFourierFeatures ( ...
                        Y_train, train_model.KY_par.gamma);
                otherwise
                    error ('Test_IOKR:InvalidArgument', ...
                        'Using features as output representation currently only "linear" and "gaussian" kernels are supported. Not %s.', ...
                        train_model.KY_par.type );
            end % switch
            
            [Psi_train, process_output] = output_feature_preprocessing_train (Y_train, ker_center);
        case 'kernel'
            [~, process_output] = output_kernel_preprocessing_train (Y_train, ...
                train_model.KY_par, ker_center);
    end % switch
    
    if (calculate_squared_norm_of_prediction)
        [KY_train, ~] = output_kernel_preprocessing_train (Y_train, ...
                    train_model.KY_par, ker_center);
    end % if

    % Scoring
    n_test = Y_C_test.getNumberOfExamples();
    scores = cell (n_test, 1);
    if (calculate_squared_norm_of_prediction)
        hh = NaN (n_test, 1);
    end % if
    for j = 1:n_test
        % TODO: Handle the case that the candidate set is empty.
        n_cand = Y_C_test.getCandidateSet (j, false, 'num');
        if (isnan (n_cand))
            scores{j} = NaN;
            
            continue;
        end % if
        
        switch train_model.KY_par.representation
            case 'feature'
                Y_Cj = Y_C_test.getCandidateSet (j, false, 'data');
                
                switch train_model.KY_par.type 
                    case 'linear'
                        % Y_Cj = Y_Cj
                    case 'gaussian'
                        Y_Cj = train_model.KY_par.rff.getRandomFourierFeatures ( ...
                            Y_Cj, train_model.KY_par.gamma);
                    otherwise
                        error ('Test_IOKR:InvalidArgument', ...
                            'Using features as output representation currently only "linear" and "gaussian" kernels are supported. Not %s.', ...
                            train_model.KY_par.type );
                end % switch      
                
                Psi_Cj = norma (Y_Cj, process_output.mean, ker_center);
                
                scores{j} = (Psi_train * B(:,j))' * Psi_Cj;
                
            case 'kernel'
                
                KY_train_Cj = output_kernel_preprocessing_test ( ...
                    Y_train, Y_C_test.getCandidateSet (j, false, 'data'), ...
                    train_model.KY_par, process_output, ker_center);

                scores{j} = B(:,j)' * KY_train_Cj;
        end % switch
        
        if (calculate_squared_norm_of_prediction)
            % <h(x_j),h(x_j)> = ||h(x_j)||^2
            hh(j) = B(:,j)' * KY_train * B(:,j);
        end % if
    end % for
    
    if (calculate_squared_norm_of_prediction)
        varargout{1} = hh;
    end % if
end % function


