function [ ] = run_IOKR_small (inputDir, outputDir)
%======================================================
% DESCRIPTION:
% Script for running MP-IOKR on a the small (MS/MS, CCS) dataset. 
%
% INPUTS:
% inputDir:     string, base-dictionary containing the input data
% result_dir:   string, base-dictionary where the output should be stored
%
%======================================================

    %--------------------------------------------------------------
    % Set up parameters
    %--------------------------------------------------------------
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (struct(), ...
        {'debug_param', 'opt_param', 'iokr_param', 'data_param', 'ky_param'});
    
    % Number of cross-validation folds for the hyper-parameter
    % optimization, e.g. finding lambda (regularization parameter).
    param.opt_param.nInnerFolds = 10;

    % Output kernel: kappa_y
    param.ky_param.representation  = 'kernel';
    param.ky_param.type            = 'gaussian';
    param.ky_param.base_kernel     = 'tanimoto';
    param.ky_param.param_selection = 'entropy';
    
    param.debug_param.verbose = true;
    
    %--------------------------------------------------------------
    % Load and prepare data
    %--------------------------------------------------------------
    % Information about MS/MS spectra ids, CCS values and molecular
    % formulas for the training compounds. The order of the comounds in
    % this table, corresponds to the order of the rows and columns in the
    % input kernel matrices. 
    cmps = readtable(fullfile(inputDir, 'compounds.csv'));
    identifier = cmps.inchikey1;
    
    % CCS values
    Z = cmps.CCS;
    % Calculate a kernel on the CCS values
    KZ = Z * Z';
    
    % Fingerprints
    fps_dir = fullfile(inputDir, 'fingerprints');
    [Y, ~] = loadFingerprints(fps_dir, fps_dir, cmps.spec_id);  % Y with shape = (n_fps, n_samples)
    [n_fps, n_samples] = size(Y);
    
    % Cross-validation indices
    cv_info = readtable(fullfile(inputDir, 'crossvalidation_folds.csv'));
    cv_params = struct('outer', struct ('type', 'fixed', 'cvInd', cv_info.cv_fold + 1), ... 
                       'inner', struct ('nFolds', param.opt_param.nInnerFolds));
    cv = getCVIndices(cv_params);
    param.opt_param.nOuterFolds = cv.outer.NumTestSets;
    n_folds = param.opt_param.nOuterFolds;
    
    
    % Molecular candidate sets
    Y_C = CandidateSetsFile_2(fullfile(inputDir, 'crossvalidation_candidates'), ...
        cmps.molecular_formula, true(n_fps, 1)); 
    assert (Y_C.getNumberOfExamples() == size (Y, 2))
         
    % Input kernels
    kernel_dir = fullfile(inputDir, 'kernels');
    kernel_files = dir(fullfile(kernel_dir, '*.txt'));
    param.data_param.availInputKernels = arrayfun(@(file) basename(file.name), ...
        kernel_files, 'UniformOutput', false);
    KX_list = loadInputKernelsIntoList(kernel_dir , param, '.txt');
    
    %--------------------------------------------------------------
    % Run Cross-validation
    %--------------------------------------------------------------
    rank = NaN(n_samples, 1);
    cand_num = arrayfun(@(id) Y_C.getCandidateSet(id, false, 'num'), ...
        1:Y_C.getNumberOfExamples());
    for i = 1:n_folds
        disp(['Now starting iteration ', int2str(i), ' out of ', int2str(n_folds)])

        % Create training and test sets
        test_set  = find(test_my(cv.outer, i));
        train_set = setdiff (1:n_samples, test_set);
        
        % Get training input and output data:
        KX_list_train = cellfun(@(x) x(train_set, train_set), KX_list, 'UniformOutput', false);
        KZ_train      = KZ(train_set, train_set);
        Y_train       = Y(:, train_set);

        train_model = Train_IOKR (KX_list_train, Y_train, ...
            param.ky_param, param.opt_param, param.iokr_param, param.debug_param.verbose);

        % Prediction and scoring
        KX_list_train_test = cellfun(@(x) x(train_set, test_set), KX_list, 'UniformOutput', false);
        KX_list_test       = cellfun(@(x) x(test_set, test_set),  KX_list, 'UniformOutput', false);
        Y_C_test           = Y_C.getSubset(test_set);

        [scores, ~] = Test_IOKR(KX_list_train_test, KX_list_test, train_model, ...
            Y_train, Y_C_test, param.iokr_param.center);

        % Computation of the ranks
        rank(test_set) = getRanksBasedOnScores(Y_C_test, identifier(test_set), scores);
        
        if (param.debug_param.verbose)
            rank_perc     = getRankPerc(rank, cand_num);
            rank_perc_100 = rank_perc(1:100);
            disp (round (rank_perc_100([1, 5, 10, 20]), 3));
        end % if
    end

    % Computation of the percentage of identified metabolites in the top k
    rank_perc     = getRankPerc(rank, cand_num);
    rank_perc_100 = rank_perc(1:100);

    disp(round(rank_perc_100([1, 5, 10, 20]), 3));

    filename = [outputDir, '/iokr/', 'rank_mkl=' param.iokr_param.mkl ...
        '_kernel=' param.ky_param.type ...
        '_base=' param.ky_param.base_kernel '_' param.ky_param.param_selection];
    save(filename,'rank_perc_100','-ascii');
end