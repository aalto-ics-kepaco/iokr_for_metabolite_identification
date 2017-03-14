function [ score ] = Test_IOKR( KX_list_train_test, KX_list_test, train_model, Y_train, Y_C_test, iokr_param )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

    n_test = length(Y_C_test); % number of test examples

    % Computation of the input kernel between training and test examples
    KX_train_test = input_kernel_preprocessing_test(KX_list_train_test, KX_list_test, train_model.process_input, iokr_param.center);
    
    % Prediction on the test set
    B  = model.C \ KX_train_test;
    
    % Output processing
    switch train_model.representation
        case 'feature'
            [Psi_train, process_output] = output_feature_preprocessing_train(Y_train, iokr_param.center);
        case 'kernel'
            [~, process_output] = output_kernel_preprocessing_train(Y_train, train_model.KY_par, iokr_param.center);
    end

    % Pre-image
    score = cell(n_test,1);
    for j = 1:n_test
        
        switch train_model.representation
            case 'feature'
                
                Psi_Cj = output_feature_preprocessing_test(Y_C_test{j}, process_output, iokr_param.center);
                
                score{j} = (Psi_train * B(:,j))' * Psi_Cj;
                
            case 'kernel'
                
                KY_train_Cj = output_kernel_preprocessing_test(Y_train, Y_C_test{j}, train_model.KY_par, process_output, iokr_param.center);

                score{j} = B(:,j)' * KY_train_Cj;
        end
    end

end

