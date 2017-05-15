function [ KX_train_test ] = input_kernel_preprocessing_test( KX_list_train_test, KX_list_test, train_process, ker_center )
            
    KX_train_test = zeros(size(KX_list_train_test{1}));
    
    for i = 1:length(KX_list_train_test)
        % Centering and normalization
        KX_train_test_i_cn = kernel_preprocessing_test(KX_list_train_test{i}, KX_list_test{i}, train_process(i), ker_center);
        
        KX_train_test = KX_train_test + train_process(i).w * KX_train_test_i_cn;
    end
    
end

