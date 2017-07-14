function IOKR_MP_feat_independent_test (inputDir, outputDir, param)
%% IOKR_MP_EVALUATION evaluation of the metabolite-identification using IOKR-MP
%    ALGORITHM:
%       = Load data associated with the metabolites = 
%           (1) load input-kernels (--> list of kernels for MKL)
%           (2) load fingerprints
%           (3) load candidate-sets
%
%       = Precalculate / load scenarios and statistics =
%           SEE ALGO

% create a candidate selection
% create cross-validation indices (outer and inner)
% pre-calculate the covariance matricies (using cv indices)

% LOOP outer cross-validation
%   

% ALGO
%   - provide selection paratemer
%   - provide cv-partition for the outer fold
%
%   - IF    pre-calculated data should be used --> load data from drive
%   - ELSE  calcuate: cv-partition, selection and statistic

    %% Check the input arguments and set defaults
    if (nargin < 2)
        error ('IOKR_MP_reverse_feat_evaluation:InvalidInput', ...
            'Not enough input arguments.');
    end % if
    if (nargin < 3)
        param = struct ();
    end % if   
    if (~ exist (inputDir, 'dir'))
        error ('IOKR_MP_reverse_feat_evaluation:InvalidInput', '%s: No such directory.', ...
            inputDir);
    end % if
    if (~ exist (outputDir, 'dir'))
        error ('IOKR_MP_reverse_feat_evaluation:InvalidInput', '%s: No such directory.', ...
            outputDir);
    end % if
    
    % Set the defaults values for the parameter in PARAM if needed.
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (param, {'debug_param', 'opt_param', 'mp_iokr_param', 'data_param'});    
    
    %% Load data 
    % ... input-kernels for the training examples
    kernel_files = dir ([inputDir '/kernels/*.txt']);
    param.data_param.availInputKernels = arrayfun (@(file) basename (file.name), ...
        kernel_files, 'UniformOutput', false);
    param.data_param.inputKernel = 'unimkl';
    KX_list = loadInputKernelsIntoList ([inputDir, '/kernels/'], param, '.txt');
    if (isempty (KX_list))
        error ('IOKR_MP_reverse_feat_evaluation:InvalidInput', ...
            'No kernel loaded.');
    end % if
    
    % Set seed in order to be able to reproduce the settings.
    % By default this is set to 'shuffle'. This means, that the current
    % time is used as seed. However, if we run the script on triton the
    % jobs might start at the same time. We therefore need to be able to
    % set the seed manualy. For example by using the slurm-job-id. This can
    % be done by the calling sbatch-file.
    rng (param.debug_param.randomSeed);
       
    
    % ... inchi keys, molecular formulas, fingerprints
    load ([inputDir '/compound_info.mat'], 'dt_inchi_mf_fp');
    
    % ... identifier (inchis) for the training examples
    inchis = dt_inchi_mf_fp.inchi_key_1; 
   
    % ... fingerprints for the training examples
    Y = full (dt_inchi_mf_fp.fp_masked)';
    [~,n] = size(Y);

    % Candidates description
   
    tic;
    if (~ isempty (param.debug_param.cand))
        warning ('Candidate taken from PARAM.DEBUG_PARAM.CAND');
        
        mf_corres = get_mf_corres (dt_inchi_mf_fp.molecular_formula, ...
             param.debug_param.cand);
        
        Y_C = CandidateSets (DataHandle (param.debug_param.cand), mf_corres);
        
        param.debug_param.cand = []; 
    else
        load ([inputDir '/cand.mat'], 'cand');
        
        mf_corres = get_mf_corres (dt_inchi_mf_fp.molecular_formula, cand);
        
        Y_C = CandidateSets (DataHandle (cand), mf_corres);
        
        clear cand;
    end 
    toc;
    
    % eval = find (arrayfun (@(x) cand(x).num < 3000, mf_corres));
    eval_set = true (numel (mf_corres), 1);
    
    assert (Y_C.getNumberOfExamples() == size (Y, 2), ...
        'The number of examples in the associated with the candidate sets is different from the number of example fingerprint vectors.');    
    
    %% Load / Store pre-calculated statistics  
    % Select a subset of candidates. By using the default values ALL the 
    % candidates are selected.
    cv_param = struct ('nObservations', n, 'outer', struct ('type', 'random', 'nFolds', param.opt_param.nOuterFolds));
    param.data_param.cv_param = cv_param;
    param.data_param.cv       = getCVIndices (param.data_param.cv_param);
    
    selec = getCandidateSelection (Y_C, inchis, param.data_param.selection_param);      
    Y_C.setSelectionsOfCandidateSets (selec);
    
    %% Evaluate the performance using 10-fold cv
    ranks      = NaN (Y_C.getNumberOfExamples(), 1);
    candNum    = arrayfun (@(idx) Y_C.getCandidateSet (idx, false, 'num'), 1:Y_C.getNumberOfExamples());
    debug_info = struct('lambda_opt', [], 'gamma_opt', [], 'mp_err', []);

    for foldIdx = 1:param.data_param.cv.outer.NumTestSets
        fprintf ('Outer fold: %d/%d\n', foldIdx, param.data_param.cv.outer.NumTestSets);
              
        data_param_fold = struct();
        data_param_fold.train_set = training_my (param.data_param.cv.outer, foldIdx);
        data_param_fold.test_set  = test_my (param.data_param.cv.outer, foldIdx);
        data_param_fold.usePreCalcStat = param.data_param.usePreCalcStat;
        
        % Calculate the scores for each candidate corresponding the current
        % test-examples
        [scores_test, debug_info(foldIdx)] = MP_IOKR_reverse_feat (KX_list, Y(:, data_param_fold.train_set), Y_C, ...
            param.opt_param, param.mp_iokr_param, data_param_fold, param.debug_param);
        assert (numel (scores_test) == numel (find (data_param_fold.test_set)), ...
            'There must be a score for each test-examples.');
        
        % Computation of the ranks of the test examples
        Y_C_test      = Y_C.getSubset (data_param_fold.test_set);
        inchis_test   = inchis(data_param_fold.test_set);
        eval_set_test = eval_set(data_param_fold.test_set);
        
        ranks_test = getRanksBasedOnScores (Y_C_test, inchis_test, scores_test, eval_set_test);
        

        
        ranks(data_param_fold.test_set) = ranks_test;
        
        if (param.debug_param.verbose)
            rank_perc     = getRankPerc (ranks, candNum);
            rank_perc_100 = rank_perc(1:100);
            disp (round (rank_perc_100([1, 5, 10, 20]), 3));
        end % if
        
        clear data_param_fold;
    end % for

    assert (all (isnan (ranks) == (~ eval_set)), ...
        'The examples without propper candidate set and examples without rank (due to the absence if a candidate set) must be equal.');
    
    %% Calculate rank percentages
    
    rankPerc = getRankPerc (ranks, candNum);
    
    %% Store results
    % FIXME: This result should be stored using the hash value for the
    % corresponding setting. 
    candNumSel = arrayfun (@(idx) Y_C.getCandidateSet (idx, 1, 'num'), 1:Y_C.getNumberOfExamples());
    
    result = struct ('ranks', ranks, 'rank_perc', rankPerc, 'cand_num', candNum, ...
        'cand_num_sel', candNumSel, 'debug_info', debug_info, 'opt_param', param.opt_param, ...
        'selection_param', param.data_param.selection_param); %#ok<NASGU>
    
    settingHash = DataHash (struct (                            ...
        'cv_param',        param.data_param.cv_param,           ...
        'selection_param', param.data_param.selection_param,    ...
        'repetition',      param.data_param.repetition,         ...
        'input_kernel',    upper(param.data_param.inputKernel), ...
        'center',          param.mp_iokr_param.center,          ...
        'rev_iokr',        param.mp_iokr_param.rev_iokr));
    save (strcat (outputDir, '/', settingHash, '.mat'), 'result', '-v7.3');
    
    disp (rankPerc([1, 5, 10, 20]));
    disp ('! Ready !');
end % function 