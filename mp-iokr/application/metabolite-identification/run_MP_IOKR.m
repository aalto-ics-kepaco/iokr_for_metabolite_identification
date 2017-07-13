function [ ] = run_MP_IOKR( inputDir, outputDir, cand )
%======================================================
% DESCRIPTION:
% Script for running MP-IOKR on a small test-dataset containing ~260
% MS2-spectra. 
%
% INPUTS:
% inputDir:     directory in which the data are contained
% result_dir:   directory in which the results will be saved
% cand:         candidate set
%
% NOTES:  
%   - Evaluate only on those examples with a candidate-set size smaller
%     than X:
%   eval = find (arrayfun (@(x) cand(x).num < 3000, mf_corres));
%
%======================================================

    % We fix the random seed to be able to compare the results directly
    % with mp-iokr.
    rng (10);

    %--------------------------------------------------------------
    % Load Data
    %--------------------------------------------------------------

    % inchi keys, molecular formulas, fingerprints
    load ([inputDir '/compound_info.mat'], 'dt_inchi_mf_fp');
    inchi = dt_inchi_mf_fp.inchi_key_1; 
    
    % Extract fingerprints
    Y = full (dt_inchi_mf_fp.fp_masked)';
    [~,n] = size(Y);

    % Candidates description
    mf_corres = get_mf_corres (dt_inchi_mf_fp.molecular_formula, cand);
   
    % Input kernels
    kernel_files = dir ([inputDir '/kernels/*.txt']);
    KX_list = arrayfun (@(file) loadKernel ([inputDir '/' file.name]), kernel_files, ...
        'UniformOutput', false);
    
    % Parameters
    iokr_param = struct('center',1,'mkl','unimkl','model_representation','only_C');
    select_param = struct( ...
        'cv_type','loocv', ...
        'lambda', [1e-5 1e-4 1e-3 1e-2 5e-2 1e-1 5e-1 1 10 100]);
%     ky_param = struct( ...
%         'type','gaussian', ...
%         'base_kernel','tanimoto',...
%         'param_selection','entropy');
    ky_param = struct ('type', 'linear', 'base_kernel', 'linear', 'param_selection', 'cv');
    output_param = struct('representation','kernel','kernel_param',ky_param);

    %--------------------------------------------------------------
    % Cross-validation
    %--------------------------------------------------------------

    rank = NaN(n,1);
    cand_num = zeros(n,1); % vector containing the number of candidates for each test example

    n_folds = 10; % number of folds
    cv = getCVIndices (struct ('nObservations', n, ...
        'outer', struct ('type', 'random', 'nFolds', n_folds)));

    for i = 1:n_folds
        disp(['Now starting iteration ', int2str(i), ' out of ', int2str(n_folds)])

        % Create training and test sets
        test_set  = find (test_my (cv.outer, i));
        train_set = setdiff(1:n,test_set);
        
        % Training
        KX_list_train = cellfun(@(x) x(train_set,train_set), KX_list, 'UniformOutput', false);
        Y_train = Y(:,train_set);

        train_model = Train_IOKR(KX_list_train, Y_train, output_param, ...
            select_param, iokr_param);

        % Prediction and scoring
        KX_list_train_test = cellfun(@(x) x(train_set,test_set), KX_list, 'UniformOutput', false);
        KX_list_test = cellfun(@(x) x(test_set,test_set), KX_list, 'UniformOutput', false);
        Y_C_test = arrayfun(@(x) x.data, cand(mf_corres(test_set)),'UniformOutput',false);

        score = Test_IOKR(KX_list_train_test, KX_list_test, train_model, ...
            Y_train, Y_C_test, iokr_param.center);

        % Computation of the ranks
        for j = 1:length(test_set)
            k = test_set(j);

            inchi_c = cand(mf_corres(k)).id;
            [~,IX] = sort(score{j},'descend');
            
            rank(k) = find(strcmp(inchi_c(IX), inchi{k}));
            cand_num(k) = length(score{j});
        end
    end

    % Computation of the percentage of identified metabolites in the top k
    rank_perc     = getRankPerc (rank, cand_num);
    rank_perc_100 = rank_perc(1:100);

    disp (round (rank_perc_100([1, 5, 10, 20]), 3));

    filename = [outputDir, 'rank_mkl=' iokr_param.mkl ...
        '_kernel=' ky_param.type ...
        '_base=' ky_param.base_kernel '_' ky_param.param_selection ...
        '_model_representation=', iokr_param.model_representation];
    save(filename,'rank_perc_100','-ascii');
    
end