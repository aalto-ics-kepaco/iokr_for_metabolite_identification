function [ ] = run_IOKR( inputDir, outputDir, cand, model_representation )
%======================================================
% DESCRIPTION:
% Script for running IOKR
%
% INPUTS:
% inputDir:     directory in which the data are contained
% result_dir:   directory in which the results will be saved
%
%======================================================

    rng (10);
%     addpath(genpath('..'));

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
    % load([inputDir '/cand.mat'],'cand');
    
    mf_corres = get_mf_corres (dt_inchi_mf_fp.molecular_formula, cand);
    
    % eval = find (arrayfun (@(x) cand(x).num < 3000, mf_corres));
    eval = 1:n;

    % Input kernels
    kernel_files = dir ([inputDir '/*.txt']);
    KX_list = arrayfun (@(file) loadKernel ([inputDir '/' file.name]), kernel_files, ...
        'UniformOutput', false);
    
    % Parameters
    iokr_param = struct('center',1,'mkl','unimkl');
    select_param = struct( ...
        'cv_type','loocv', ...
        'lambda', [1e-5 1e-4 1e-3 1e-2 5e-2 1e-1 5e-1 1 10 100]);
    ky_param = struct( ...
        'type','gaussian', ...
        'base_kernel','tanimoto',...
        'param_selection','entropy');
    output_param = struct('representation','kernel','kernel_param',ky_param);

    %--------------------------------------------------------------
    % Cross-validation
    %--------------------------------------------------------------

    rank = NaN(n,1);
    cand_num = zeros(n,1); % vector containing the number of candidates for each test example

    n_folds = 10; % number of folds
    % ind_fold = load([inputDir 'cv_ind.txt']); % indices of the different folds
    cv = cvpartition (n, 'Kfold', n_folds);

    for i = 1:n_folds
        disp(['Now starting iteration ', int2str(i), ' out of ', int2str(n_folds)])

        % Create training and test sets
        test_set = find (test (cv, i));
        % test_set = find(ind_fold == i);
        train_set = setdiff(1:n,test_set);
        test_set = intersect(test_set,eval);
        
        % Training
        KX_list_train = cellfun(@(x) x(train_set,train_set), KX_list, 'UniformOutput', false);
        Y_train = Y(:,train_set);

        disp ('TRAINING')
        train_model = Train_IOKR(KX_list_train, Y_train, output_param, ...
            select_param, iokr_param, model_representation);

        % Prediction 
        KX_list_train_test = cellfun(@(x) x(train_set,test_set), KX_list, 'UniformOutput', false);
        KX_list_test = cellfun(@(x) x(test_set,test_set), KX_list, 'UniformOutput', false);
        Y_C_test = arrayfun(@(x) x.data', cand(mf_corres(test_set)),'UniformOutput',false);

        disp ('TESTING')
        score = Test_IOKR(KX_list_train_test, KX_list_test, train_model, ...
            Y_train, Y_C_test, iokr_param.center, model_representation);

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
    nel = hist(rank(eval), 1:max(cand_num));
    rank_perc = cumsum(nel)';
    rank_perc = rank_perc/length(eval)*100;
    rank_perc_100 = rank_perc(1:100);
    
    disp (round (rank_perc_100([1, 5, 10, 20]), 3));

    filename = [outputDir 'rank_mkl=' iokr_param.mkl ...
        '_kernel=' ky_param.type ...
        '_base=' ky_param.base_kernel '_' ky_param.param_selection ...
        '_model_representation=', model_representation];
    save(filename,'rank_perc_100','-ascii');
    
end


