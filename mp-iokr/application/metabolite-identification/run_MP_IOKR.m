function [ ] = run_MP_IOKR (inputDir, outputDir, cand)
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

    %--------------------------------------------------------------
    % Set up parameters
    %--------------------------------------------------------------
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (struct(), ...
        {'debug_param', 'opt_param', 'mp_iokr_param', 'data_param', 'ky_param'});
    
    param.debug_param.randomSeed = 10;
    rng (param.debug_param.randomSeed);
    
    n_folds = param.opt_param.nOuterFolds;
    param.opt_param.nInnerFolds = 10;
   
    %--------------------------------------------------------------
    % Load and prepare data
    %--------------------------------------------------------------

    % inchi keys, molecular formulas, fingerprints
    load ([inputDir '/compound_info.mat'], 'dt_inchi_mf_fp');
    inchi = dt_inchi_mf_fp.inchi_key_1; 
    
    % Extract fingerprints
    Y = full (dt_inchi_mf_fp.fp_masked)';
    [~,n] = size(Y);
    param.ky_param.representation  = 'kernel';
    param.ky_param.type            = 'gaussian';
    param.ky_param.base_kernel     = 'linear';
    param.ky_param.param_selection = 'entropy';
    param.ky_param.rff_dimension   = 1000;
    

    % Cross-validation
    cv = getCVIndices (struct ('nObservations', n, ...
        'outer', struct ('type', 'random', 'nFolds', n_folds)));

    % Candidates description
    mf_corres = get_mf_corres (dt_inchi_mf_fp.molecular_formula, cand);
    Y_C       = CandidateSets (DataHandle (cand), mf_corres);
    assert (Y_C.getNumberOfExamples() == size (Y, 2))
    % Candidate selection
    param.data_param.selection_param = struct ( ...
        'strategy', 'random', 'perc', 1, 'inclExpCand', false);
    selec = getCandidateSelection (Y_C, inchi, ...
       param.data_param.selection_param);      
    Y_C.setSelectionsOfCandidateSets (selec);
    
    % Input kernels
    kernel_files = dir ([inputDir '/kernels/*.txt']);
    param.data_param.availInputKernels = arrayfun (@(file) basename (file.name), ...
        kernel_files, 'UniformOutput', false);
    KX_list = loadInputKernelsIntoList ([inputDir, '/kernels/'], param, '.txt');
    
    %--------------------------------------------------------------
    % Run Cross-validation
    %--------------------------------------------------------------
    rank     = NaN (n, 1);
    cand_num = arrayfun (@(id) Y_C.getCandidateSet (id, false, 'num'), ...
        1:Y_C.getNumberOfExamples());
    for i = 1:n_folds
        disp(['Now starting iteration ', int2str(i), ' out of ', int2str(n_folds)])

        % Create training and test sets
        test_set  = find (test_my (cv.outer, i));
        train_set = setdiff (1:n, test_set);
        
        % Training
        KX_list_train = cellfun(@(x) x(train_set, train_set), KX_list, 'UniformOutput', false);
        Y_train       = Y(:, train_set);
        Y_C_train     = Y_C.getSubset (train_set);

        train_model = Train_MPIOKR (KX_list_train, Y_train, Y_C_train, ...
            param.ky_param, param.mp_iokr_param, param.opt_param, param.debug_param);

        % Prediction and scoring
        KX_list_train_test = cellfun(@(x) x(train_set,test_set), KX_list, 'UniformOutput', false);
        KX_list_test       = cellfun(@(x) x(test_set,test_set),  KX_list, 'UniformOutput', false);
        Y_C_test           = Y_C.getSubset (test_set);

        scores = Test_MPIOKR (KX_list_train_test, KX_list_test, train_model, ...
            Y_train, Y_C_train, Y_C_test, param.mp_iokr_param, param.mp_iokr_param.center, param.debug_param);

        % Computation of the ranks
        rank(test_set) = getRanksBasedOnScores (Y_C_test, inchi(test_set), scores);
        
        if (param.debug_param.verbose)
            rank_perc     = getRankPerc (rank, cand_num);
            rank_perc_100 = rank_perc(1:100);
            disp (round (rank_perc_100([1, 5, 10, 20]), 3));
        end % if
    end

    % Computation of the percentage of identified metabolites in the top k
    rank_perc     = getRankPerc (rank, cand_num);
    rank_perc_100 = rank_perc(1:100);

    disp (round (rank_perc_100([1, 5, 10, 20]), 3));

    filename = [outputDir, '/mpiokr/', 'rank_mkl=', param.mp_iokr_param.mkl, ...
        '_kernel=', param.ky_param.type, ...
        '_base=', param.ky_param.base_kernel, ...
        '_', param.ky_param.param_selection, ...
        '_strategy=', param.data_param.selection_param.strategy, ...
        '_inclCandExp', num2str(param.data_param.selection_param.inclExpCand)];
    save(filename,'rank_perc_100','-ascii');
end