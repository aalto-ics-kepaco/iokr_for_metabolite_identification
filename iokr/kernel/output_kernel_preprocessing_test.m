function [ KY_train_cand_cn ] = output_kernel_preprocessing_test( Y_train, Y_cand, KY_par, train_process, ker_center )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

    KY_train_cand = build_kernel(Y_train, Y_cand, KY_par);
    KY_cand = build_kernel(Y_cand, Y_cand, KY_par);
    
    KY_train_cand_cn = kernel_preprocessing_test(KY_train_cand, KY_cand, train_process, ker_center);


end

