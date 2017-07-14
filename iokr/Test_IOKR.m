function [ score, process_output ] = Test_IOKR( KX_list_train_test, KX_list_test, ...
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

    % Computation of the input kernel between training and test examples
    KX_train_test = input_kernel_preprocessing_test(KX_list_train_test, ...
        KX_list_test, train_model.process_input, ker_center);
    
    % Prediction on the test set 
    switch train_model.model_representation
        case 'only_C'
            B = train_model.C \ KX_train_test;
        case 'Chol_decomp_of_C'
            y = linsolve (train_model.C,  KX_train_test, struct ('LT', true));
            B = linsolve (train_model.C', y,             struct ('UT', true));
    end % switch
    % Pre-image
    
    % Preprocessing of the training outputs    
    switch train_model.representation
        case 'feature'
            [Psi_train, process_output] = ...
                output_feature_preprocessing_train(Y_train, ker_center);
        case 'kernel'
            [~,         process_output] = ...
                output_kernel_preprocessing_train(Y_train, ...
                train_model.KY_par, ker_center);
    end
   
    % Scoring
    n_test = Y_C_test.getNumberOfExamples();
    score = cell(n_test,1);
    for j = 1:n_test    
        switch train_model.representation
            case 'feature'
                                
                Psi_Cj = norma(Y_C_test.getCandidateSet (j, false, 'data'), ...
                    process_output.mean, ker_center);
                
                score{j} = (Psi_train * B(:,j))' * Psi_Cj;
                
            case 'kernel'
                
                KY_train_Cj = output_kernel_preprocessing_test( ...
                    Y_train, Y_C_test.getCandidateSet (j, false, 'data'), ...
                    train_model.KY_par, process_output, ker_center);

                score{j} = B(:,j)' * KY_train_Cj;
        end
    end
end


